"""Espelho independente Python/SciPy dos estudos paramétricos do OpenOtto v0.3."""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


PI = math.pi
RU = 8.31446261815324
T_REF = 298.15
P_INT = 101.35e3
T_INT = 303.15
P_EXH = 110.0e3
T_EXH = 750.0
T_WALL = 430.0
H_WALL = 350.0
FRICTION_COEFFICIENT = 16.0

VD_TOTAL = 1667e-6
CYLINDERS = 6
VD = VD_TOTAL / CYLINDERS
BORE = (4.0 * VD / PI) ** (1.0 / 3.0)
STROKE = BORE
R_CRANK = STROKE / 2.0
ROD_LENGTH = 3.7 * R_CRANK
PISTON_AREA = PI * BORE**2 / 4.0

MM = np.array([16.04246, 31.9988, 28.0134, 44.0095, 18.01528, 28.0101, 2.01588]) * 1e-3
NU = np.array([-1.0, -2.0, 0.0, 1.0, 2.0, 0.0, 0.0])

THERMO = [
    ((5.14987613, -1.36709788e-2, 4.91800599e-5, -4.84743026e-8, 1.66693956e-11, -1.02466476e4, -4.64130376), (7.48514950e-2, 1.33909467e-2, -5.73285809e-6, 1.22292535e-9, -1.01815230e-13, -9.46834459e3, 18.4373180)),
    ((3.78245636, -2.99673416e-3, 9.84730201e-6, -9.68129509e-9, 3.24372837e-12, -1.06394356e3, 3.65767573), (3.28253784, 1.48308754e-3, -7.57966669e-7, 2.09470555e-10, -2.16717794e-14, -1.08845772e3, 5.45323129)),
    ((3.53100528, -1.23660987e-4, -5.02999433e-7, 2.43530612e-9, -1.40881235e-12, -1.04697628e3, 2.96747468), (2.95257626, 1.39690040e-3, -4.92631603e-7, 7.86010195e-11, -4.60755204e-15, -9.23948645e2, 5.87188762)),
    ((2.35677352, 8.98459677e-3, -7.12356269e-6, 2.45919022e-9, -1.43699548e-13, -4.83719697e4, 9.90105222), (3.85796028, 4.41437026e-3, -2.21481404e-6, 5.23490188e-10, -4.72084164e-14, -4.87591660e4, 2.27163806)),
    ((4.19864056, -2.03643410e-3, 6.52040211e-6, -5.48797062e-9, 1.77197817e-12, -3.02937267e4, -0.849032208), (3.03399249, 2.17691804e-3, -1.64072518e-7, -9.70419870e-11, 1.68200992e-14, -3.00042971e4, 4.96677010)),
    ((3.57953347, -6.10353680e-4, 1.01681433e-6, 9.07005884e-10, -9.04424499e-13, -1.43440860e4, 3.50840928), (2.71518561, 2.06252743e-3, -9.98825771e-7, 2.30053008e-10, -2.03647716e-14, -1.41518724e4, 7.81868772)),
    ((2.34433112, 7.98052075e-3, -1.94781510e-5, 2.01572094e-8, -7.37611761e-12, -9.17935173e2, 0.683010238), (3.33727920, -4.94024731e-5, 4.99456778e-7, -1.79566394e-10, 2.00255376e-14, -9.50158922e2, -3.20502331)),
]


@dataclass(frozen=True)
class ModelConfig:
    rpm: float = 2000.0
    phi: float = 1.0
    compression_ratio: float = 13.0
    property_model: str = "nasa"


def piston_position(angle: np.ndarray | float) -> np.ndarray | float:
    return R_CRANK * (1.0 - np.cos(angle)) + ROD_LENGTH - np.sqrt(
        ROD_LENGTH**2 - R_CRANK**2 * np.sin(angle) ** 2
    )


def volume(angle: np.ndarray | float, cfg: ModelConfig) -> np.ndarray | float:
    clearance = VD / (cfg.compression_ratio - 1.0)
    return clearance + PISTON_AREA * piston_position(angle)


def dvolume(angle: np.ndarray | float) -> np.ndarray | float:
    sine = np.sin(angle)
    cosine = np.cos(angle)
    return PISTON_AREA * (
        R_CRANK * sine
        + R_CRANK**2 * sine * cosine
        / np.sqrt(ROD_LENGTH**2 - R_CRANK**2 * sine**2)
    )


def chamber_area(angle: np.ndarray | float) -> np.ndarray | float:
    return 2.0 * PISTON_AREA + PI * BORE * piston_position(angle)


def valve_area(angle: float, opens: float, closes: float, maximum: float) -> float:
    degrees = math.degrees(angle)
    for shift in (-1, 0, 1):
        shifted = degrees + 720.0 * shift
        if opens < shifted < closes:
            phase = (shifted - opens) / (closes - opens)
            return maximum * math.sin(PI * phase) ** 2
    return 0.0


def wiebe_derivative(angle: float) -> float:
    z = (math.degrees(angle) + 15.0) / 50.0
    if z <= 0.0 or z >= 1.0:
        return 0.0
    return 15.0 * z**2 * math.exp(-5.0 * z**3) / (
        math.radians(50.0) * (1.0 - math.exp(-5.0))
    )


def nasa_coeff(index: int, temperature: float):
    return THERMO[index][0 if temperature <= 1000.0 else 1]


def nasa_cp(index: int, temperature: float) -> float:
    a = nasa_coeff(index, temperature)
    return RU * (
        a[0] + a[1] * temperature + a[2] * temperature**2
        + a[3] * temperature**3 + a[4] * temperature**4
    )


def nasa_h(index: int, temperature: float) -> float:
    a = nasa_coeff(index, temperature)
    return RU * temperature * (
        a[0] + a[1] * temperature / 2.0 + a[2] * temperature**2 / 3.0
        + a[3] * temperature**3 / 4.0 + a[4] * temperature**4 / 5.0
        + a[5] / temperature
    )


def species_cp(cfg: ModelConfig, index: int, temperature: float) -> float:
    return nasa_cp(index, temperature if cfg.property_model == "nasa" else T_REF)


def species_h(cfg: ModelConfig, index: int, temperature: float) -> float:
    if cfg.property_model == "nasa":
        return nasa_h(index, temperature)
    return nasa_h(index, T_REF) + nasa_cp(index, T_REF) * (temperature - T_REF)


def mixture_u(cfg: ModelConfig, moles: np.ndarray, temperature: float) -> float:
    return sum(
        max(moles[i], 0.0) * (species_h(cfg, i, temperature) - RU * temperature)
        for i in range(7)
    )


def temperature_from_u(cfg: ModelConfig, moles: np.ndarray, energy: float) -> float:
    lower, upper = 150.0, 6500.0
    for _ in range(42):
        middle = (lower + upper) / 2.0
        if mixture_u(cfg, moles, middle) > energy:
            upper = middle
        else:
            lower = middle
    return (lower + upper) / 2.0


def intake_composition(cfg: ModelConfig) -> np.ndarray:
    amounts = np.array([1.0, 2.0 / cfg.phi, 7.52 / cfg.phi, 0.0, 0.0, 0.0, 0.0])
    return amounts / amounts.sum()


def exhaust_composition(cfg: ModelConfig) -> np.ndarray:
    amounts = np.array([0.0, 2.0 / cfg.phi - 2.0, 7.52 / cfg.phi, 1.0, 2.0, 0.0, 0.0])
    return amounts / amounts.sum()


def orifice_mass_flow(area, cd, p_up, t_up, p_down, r_specific, gamma):
    if area <= 0.0 or p_up <= p_down:
        return 0.0
    ratio = np.clip(p_down / p_up, 0.0, 1.0)
    critical = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    if ratio <= critical:
        factor = math.sqrt(
            gamma * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (gamma - 1.0))
        )
    else:
        term = 2.0 * gamma / (gamma - 1.0) * (
            ratio ** (2.0 / gamma) - ratio ** ((gamma + 1.0) / gamma)
        )
        factor = math.sqrt(max(term, 0.0))
    return cd * area * p_up / math.sqrt(r_specific * t_up) * factor


def port_flow(cfg, angle, p_cyl, t_cyl, moles, reservoir, valve, x_res):
    p_res, t_res = reservoir
    opens, closes, maximum, cd = valve
    area = valve_area(angle, opens, closes, maximum)
    x_cyl = np.maximum(moles, 0.0)
    x_cyl /= x_cyl.sum()
    if p_res >= p_cyl:
        x_up, t_up, p_up, p_down, direction = x_res, t_res, p_res, p_cyl, 1.0
    else:
        x_up, t_up, p_up, p_down, direction = x_cyl, t_cyl, p_cyl, p_res, -1.0
    cp_mix = sum(x_up[i] * species_cp(cfg, i, t_up) for i in range(7))
    gamma = cp_mix / (cp_mix - RU)
    molar_mass = float(x_up @ MM)
    mass_rate = orifice_mass_flow(
        area, cd, p_up, t_up, p_down, RU / molar_mass, gamma
    )
    molar_rate = mass_rate / molar_mass
    species_rate = direction * molar_rate * x_up
    enthalpy_rate = direction * molar_rate * sum(
        x_up[i] * species_h(cfg, i, t_up) for i in range(7)
    )
    return species_rate, direction * mass_rate, enthalpy_rate


def reaction_heat() -> float:
    return -sum(NU[i] * nasa_h(i, T_REF) for i in range(7))


Q_MOL = reaction_heat()


def rhs(cfg: ModelConfig):
    omega = 2.0 * PI * cfg.rpm / 60.0
    x_int = intake_composition(cfg)
    x_exh = exhaust_composition(cfg)

    def evaluate(angle: float, state: np.ndarray) -> np.ndarray:
        moles = state[:7]
        temperature = temperature_from_u(cfg, moles, state[7])
        pressure = np.maximum(moles, 0.0).sum() * RU * temperature / volume(angle, cfg)
        species_i, mdot_i, hdot_i = port_flow(
            cfg, angle, pressure, temperature, moles,
            (P_INT, T_INT), (-380.0, -130.0, 3.0e-4, 0.70), x_int,
        )
        species_e, mdot_e, hdot_e = port_flow(
            cfg, angle, pressure, temperature, moles,
            (P_EXH, T_EXH), (130.0, 380.0, 2.5e-4, 0.72), x_exh,
        )
        result = np.zeros(16)
        result[:7] = (species_i + species_e) / omega
        burn_rate = wiebe_derivative(angle)
        if burn_rate > 0.0 and state[15] > 0.0:
            extent_rate = state[15] * burn_rate
            result[:7] += NU * extent_rate
            result[10] = Q_MOL * extent_rate
        dv = dvolume(angle)
        wall_heat = H_WALL * chamber_area(angle) * (temperature - T_WALL) / omega
        indicated_work = pressure * dv
        result[7] = -wall_heat - indicated_work + (hdot_i + hdot_e) / omega
        result[8] = indicated_work
        result[9] = FRICTION_COEFFICIENT * omega * (dv / PISTON_AREA) ** 2
        result[11] = wall_heat
        result[12] = (hdot_i + hdot_e) / omega
        result[13] = max(mdot_i, 0.0) / omega
        result[14] = max(-mdot_e, 0.0) / omega
        return result

    return evaluate


def integrate_cycle(cfg: ModelConfig, base: np.ndarray, t_eval=None):
    ignition = math.radians(-15.0)
    state0 = np.zeros(16)
    state0[:8] = base
    derivative = rhs(cfg)
    atol = np.r_[np.ones(7) * 1e-11, 1e-7, np.ones(7) * 1e-8, 1e-11]
    eval_first = None if t_eval is None else t_eval[t_eval <= ignition]
    first = solve_ivp(
        derivative, (-2.0 * PI, ignition), state0, t_eval=eval_first,
        rtol=3e-8, atol=atol, max_step=math.radians(0.5),
    )
    at_ignition = first.y[:, -1].copy()
    at_ignition[15] = min(
        max(at_ignition[0], 0.0), max(at_ignition[1], 0.0) / 2.0
    )
    eval_second = None if t_eval is None else t_eval[t_eval > ignition]
    second = solve_ivp(
        derivative, (ignition, 2.0 * PI), at_ignition, t_eval=eval_second,
        rtol=3e-8, atol=atol, max_step=math.radians(0.5),
    )
    if t_eval is None:
        return second.y[:8, -1], None, None
    times = np.concatenate([first.t, second.t])
    states = np.concatenate([first.y, second.y], axis=1)
    return states[:8, -1], times, states


def solve_case(cfg: ModelConfig) -> dict:
    x_exh = exhaust_composition(cfg)
    total_moles = P_EXH * volume(-2.0 * PI, cfg) / (RU * T_EXH)
    moles = total_moles * x_exh
    base = np.r_[moles, mixture_u(cfg, moles, T_EXH)]
    cycles = 0
    change = math.inf
    for cycles in range(1, 61):
        new_base, _, _ = integrate_cycle(cfg, base)
        old_temperature = temperature_from_u(cfg, base[:7], base[7])
        new_temperature = temperature_from_u(cfg, new_base[:7], new_base[7])
        composition_change = np.max(np.abs(new_base[:7] - base[:7])) / max(
            new_base[:7].sum(), 1e-12
        )
        old_mass = float(np.maximum(base[:7], 0.0) @ MM)
        new_mass = float(np.maximum(new_base[:7], 0.0) @ MM)
        change = max(
            composition_change,
            abs(new_mass - old_mass) / max(new_mass, 1e-12),
            abs(new_temperature - old_temperature) / max(new_temperature, 1e-12),
        )
        base = new_base
        if cycles >= 3 and change <= 2e-6:
            break

    grid = np.radians(np.arange(-360.0, 360.0001, 1.0))
    cycle_initial = base.copy()
    final, times, states = integrate_cycle(cfg, cycle_initial, grid)
    temperatures = np.array([
        temperature_from_u(cfg, states[:7, i], states[7, i])
        for i in range(states.shape[1])
    ])
    pressures = (
        np.maximum(states[:7], 0.0).sum(axis=0) * RU * temperatures
        / volume(times, cfg)
    )
    indicated_work = states[8, -1]
    friction_work = states[9, -1]
    released_energy = states[10, -1]
    wall_heat = states[11, -1]
    enthalpy_flow = states[12, -1]
    delta_u = final[7] - cycle_initial[7]
    energy_residual = delta_u - (-wall_heat - indicated_work + enthalpy_flow)
    reference_molar_mass = float(intake_composition(cfg) @ MM)
    reference_mass = P_INT * reference_molar_mass * VD / (RU * T_INT)

    return {
        "config": asdict(cfg),
        "cycles": cycles,
        "periodic_change": change,
        "alpha_deg": np.degrees(times),
        "volume_cm3": volume(times, cfg) * 1e6,
        "pressure_MPa": pressures * 1e-6,
        "temperature_K": temperatures,
        "summary": {
            "released_energy_J": released_energy,
            "indicated_work_J": indicated_work,
            "friction_work_J": friction_work,
            "brake_work_J": indicated_work - friction_work,
            "eta_i": indicated_work / released_energy,
            "eta_b": (indicated_work - friction_work) / released_energy,
            "imep_kPa": indicated_work / VD / 1e3,
            "bmep_kPa": (indicated_work - friction_work) / VD / 1e3,
            "volumetric_efficiency": states[13, -1] / reference_mass,
            "pmax_MPa": float(np.max(pressures) * 1e-6),
            "Tmax_K": float(np.max(temperatures)),
            "wall_heat_J": wall_heat,
            "energy_residual_J": energy_residual,
        },
    }


STUDIES = {
    "rotation": {
        "title": "Efeito da rotação",
        "cases": [
            ("2000 rpm", ModelConfig(rpm=2000.0)),
            ("4000 rpm", ModelConfig(rpm=4000.0)),
            ("6000 rpm", ModelConfig(rpm=6000.0)),
            ("8000 rpm", ModelConfig(rpm=8000.0)),
        ],
    },
    "mixture": {
        "title": "Efeito da composição da mistura",
        "cases": [
            ("φ = 0,70", ModelConfig(phi=0.70)),
            ("φ = 0,85", ModelConfig(phi=0.85)),
            ("φ = 1,00", ModelConfig(phi=1.00)),
        ],
    },
    "compression": {
        "title": "Efeito da razão de compressão",
        "cases": [
            ("r = 8", ModelConfig(compression_ratio=8.0)),
            ("r = 10", ModelConfig(compression_ratio=10.0)),
            ("r = 13", ModelConfig(compression_ratio=13.0)),
            ("r = 16", ModelConfig(compression_ratio=16.0)),
        ],
    },
    "properties": {
        "title": "Comparação das propriedades termodinâmicas",
        "cases": [
            ("cp constante", ModelConfig(property_model="constant")),
            ("cp(T), NASA", ModelConfig(property_model="nasa")),
        ],
    },
}


def plot_study(slug: str, definition: dict, results: dict, output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    colors = ["#2563a6", "#c84d32", "#2f855a", "#7c3aed"]
    linestyles = ["-", "--", "-.", ":"]
    for index, (label, cfg) in enumerate(definition["cases"]):
        result = results[cfg]
        style = {
            "color": colors[index],
            "linestyle": linestyles[index],
            "linewidth": 2.0,
            "label": label,
        }
        axes[0].plot(result["volume_cm3"], result["pressure_MPa"], **style)
        mask = (result["alpha_deg"] >= -180.0) & (result["alpha_deg"] <= 180.0)
        axes[1].plot(result["alpha_deg"][mask], result["pressure_MPa"][mask], **style)
        axes[2].plot(result["volume_cm3"], result["temperature_K"], **style)

    axes[0].set_title("P × V")
    axes[0].set_xlabel("Volume V [cm³]")
    axes[0].set_ylabel("Pressão p [MPa]")
    axes[1].set_title("P × α")
    axes[1].set_xlabel("Ângulo α [graus]")
    axes[1].set_ylabel("Pressão p [MPa]")
    axes[1].set_xlim(-180.0, 180.0)
    axes[2].set_title("T × V")
    axes[2].set_xlabel("Volume V [cm³]")
    axes[2].set_ylabel("Temperatura T [K]")
    for axis in axes:
        axis.grid(True, color="#d9dde3", linewidth=0.7)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(labels), frameon=False,
               bbox_to_anchor=(0.5, 0.91))
    fig.suptitle(definition["title"], fontsize=16, fontweight="semibold", y=0.99)
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.14, top=0.78, wspace=0.28)
    fig.savefig(output_dir / f"{slug}_pv_palpha_tv.png", dpi=240, facecolor="white")
    plt.close(fig)


def serializable_result(result: dict) -> dict:
    converted = dict(result)
    for key in ("alpha_deg", "volume_cm3", "pressure_MPa", "temperature_K"):
        converted[key] = np.round(result[key], 7).tolist()
    converted["summary"] = {key: float(value) for key, value in result["summary"].items()}
    return converted


def main() -> None:
    output_dir = Path(__file__).resolve().parents[1] / "output" / "estudos_parametricos"
    output_dir.mkdir(parents=True, exist_ok=True)
    unique_configs = []
    for definition in STUDIES.values():
        for _, cfg in definition["cases"]:
            if cfg not in unique_configs:
                unique_configs.append(cfg)

    results = {}
    for cfg in unique_configs:
        print(f"Executando {cfg}", flush=True)
        results[cfg] = solve_case(cfg)
        summary = results[cfg]["summary"]
        print(
            f"  ciclos={results[cfg]['cycles']}; pmax={summary['pmax_MPa']:.4f} MPa; "
            f"Tmax={summary['Tmax_K']:.1f} K; resíduo={summary['energy_residual_J']:.3e} J",
            flush=True,
        )

    payload = {"studies": {}}
    for slug, definition in STUDIES.items():
        plot_study(slug, definition, results, output_dir)
        payload["studies"][slug] = {
            "title": definition["title"],
            "cases": [
                {"label": label, **serializable_result(results[cfg])}
                for label, cfg in definition["cases"]
            ],
        }

    json_path = output_dir / "parametric_results.json"
    json_path.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
                         encoding="utf-8")
    print(json_path)


if __name__ == "__main__":
    main()
