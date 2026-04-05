import streamlit as st
import pandas as pd
import numpy as np
import math
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from io import BytesIO
import json

# ================== 常數與預設參數 ==================
G = 9.81
RHO = 1.2
FR = 0.015
CD = 0.4
ETA_DRIVE = 0.9                # 傳動效率
ETA_MOTOR = 1.0                 # 馬達效率（不再使用，改為使用者輸入）
ETA_CONTROLLER = 0.95
BATTERY_ENERGY_DENSITY = 150
MOTOR_POWER_DENSITY = 1.0
CELL_VOLTAGE = 3.7
CELL_CAPACITY = 2.5

# ================== 輔助函數 ==================
def get_cd_by_vehicle(vehicle_type):
    mapping = {
        '小型電動車': 0.3,
        '電動機車': 0.5,
        '電動三輪車': 0.45,
        '高爾夫球車': 0.4
    }
    return mapping.get(vehicle_type, CD)

def calculate_power_requirements(mass, speed_ms, area, cd, fr):
    F_roll = mass * G * fr
    F_air = 0.5 * RHO * cd * area * speed_ms**2
    F_total = F_roll + F_air
    P_wheel = F_total * speed_ms
    P_motor = P_wheel / ETA_DRIVE / 1000
    return P_motor, F_total

def estimate_motor_from_power(max_power_kw, voltage, n_max, motor_eff_percent, base_speed=3000):
    rated_power = max_power_kw / 2
    T_peak = (max_power_kw * 1000) / (2 * math.pi * base_speed / 60)
    T_rated = T_peak / 2
    motor_spec = {
        '類型': '永磁同步馬達 (PMSM)',
        '最大功率 (kW)': round(max_power_kw, 2),
        '額定功率 (kW)': round(rated_power, 2),
        '峰值扭矩 (Nm)': round(T_peak, 1),
        '額定扭矩 (Nm)': round(T_rated, 1),
        '額定轉速 (rpm)': base_speed,
        '最高轉速 (rpm)': round(n_max, 0),
        '基速 (rpm)': base_speed,
        '電壓 (V)': voltage,
        '效率 (%)': round(motor_eff_percent, 1),
        '估計重量 (kg)': round(max_power_kw * MOTOR_POWER_DENSITY, 1)
    }
    return motor_spec, base_speed, T_peak

def estimate_motor_from_params(max_power_kw, peak_torque_Nm, voltage, n_max, motor_eff_percent):
    base_speed = (max_power_kw * 1000 * 60) / (2 * math.pi * peak_torque_Nm)
    if base_speed > n_max:
        base_speed = n_max
    rated_power = max_power_kw / 2
    T_rated = peak_torque_Nm / 2
    motor_spec = {
        '類型': '永磁同步馬達 (PMSM)',
        '最大功率 (kW)': round(max_power_kw, 2),
        '額定功率 (kW)': round(rated_power, 2),
        '峰值扭矩 (Nm)': round(peak_torque_Nm, 1),
        '額定扭矩 (Nm)': round(T_rated, 1),
        '額定轉速 (rpm)': round(base_speed, 0),
        '最高轉速 (rpm)': round(n_max, 0),
        '基速 (rpm)': round(base_speed, 0),
        '電壓 (V)': voltage,
        '效率 (%)': round(motor_eff_percent, 1),
        '估計重量 (kg)': round(max_power_kw * MOTOR_POWER_DENSITY, 1)
    }
    return motor_spec, base_speed, peak_torque_Nm

def estimate_battery(rated_power_kw, voltage, duration_h=1.0):
    energy_kwh = rated_power_kw * duration_h
    capacity_ah = (energy_kwh * 1000) / voltage
    c_rate = 1.0 / duration_h
    series = round(voltage / CELL_VOLTAGE)
    parallel = math.ceil(capacity_ah / CELL_CAPACITY)
    weight = (energy_kwh * 1000) / BATTERY_ENERGY_DENSITY
    battery_spec = {
        '類型': '鋰離子電池 (Li-ion)',
        '標稱電壓 (V)': voltage,
        '容量 (Ah)': round(capacity_ah, 1),
        '能量 (kWh)': round(energy_kwh, 2),
        '放電倍率 (C)': round(c_rate, 1),
        '串聯數': series,
        '並聯數': parallel,
        '估計重量 (kg)': round(weight, 1)
    }
    return battery_spec

def estimate_controller(max_power_kw, voltage):
    I_max = (max_power_kw * 1000) / voltage
    controller_spec = {
        '類型': 'MOSFET 控制器',
        '最大功率 (kW)': round(max_power_kw, 2),
        '電壓範圍 (V)': f"{int(voltage*0.8)}-{int(voltage*1.2)}",
        '最大電流 (A)': round(I_max, 1),
        '效率 (%)': ETA_CONTROLLER * 100
    }
    return controller_spec

def estimate_gearbox(speed_max_ms, wheel_radius_m):
    wheel_rpm = speed_max_ms * 60 / (2 * math.pi * wheel_radius_m)
    gear_ratio = 6000 / wheel_rpm
    return gear_ratio

def calculate_load_curve(mass, area, cd, fr, wheel_radius_m, gear_ratio, speed_max_ms, grade_percent=0):
    speeds = np.linspace(0, speed_max_ms * 1.1, 100)
    torques = []
    forces = []
    grade_rad = math.atan(grade_percent / 100)
    for v in speeds:
        if v == 0:
            torques.append(0)
            forces.append(0)
            continue
        F_roll = mass * G * fr * math.cos(grade_rad)
        F_air = 0.5 * RHO * cd * area * v**2
        F_grade = mass * G * math.sin(grade_rad)
        F_total = F_roll + F_air + F_grade
        T_wheel = F_total * wheel_radius_m
        T_motor = T_wheel / gear_ratio / ETA_DRIVE
        torques.append(T_motor)
        forces.append(F_total)
    motor_rpm = speeds * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio
    return motor_rpm, np.array(torques), speeds * 3.6, np.array(forces)

def find_intersection(x1, y1, x2, y2):
    y2_interp = np.interp(x1, x2, y2)
    diff = y1 - y2_interp
    intersections = []
    for i in range(len(x1)-1):
        if diff[i] * diff[i+1] <= 0:
            x_cross = x1[i] - diff[i] * (x1[i+1] - x1[i]) / (diff[i+1] - diff[i])
            y_cross = np.interp(x_cross, x1, y1)
            intersections.append((x_cross, y_cross))
    return intersections

def simulate_acceleration(mass, area, cd, fr, wheel_radius_m, gear_ratio, motor_spec, base_speed, T_peak, speed_max_ms, dt=0.1):
    n_max = motor_spec['最高轉速 (rpm)']
    P_peak = motor_spec['最大功率 (kW)']

    def get_max_torque(v):
        if v <= 0:
            return T_peak
        n = v * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio
        if n <= base_speed:
            return T_peak
        elif n <= n_max:
            return (P_peak * 1000) / (2 * math.pi * n / 60)
        else:
            return 0

    t = 0
    v = 0
    x = 0
    time_list = [0]
    speed_list = [0]
    disp_list = [0]

    def resistance(v_ms):
        F_roll = mass * G * fr
        F_air = 0.5 * RHO * cd * area * v_ms**2
        return F_roll + F_air

    while v < speed_max_ms * 0.99 and t < 60:
        T_motor = get_max_torque(v)
        F_drive = T_motor * gear_ratio * ETA_DRIVE / wheel_radius_m
        F_resist = resistance(v)
        F_net = F_drive - F_resist
        a = F_net / mass
        if a < 0:
            break
        v += a * dt
        x += v * dt
        t += dt
        time_list.append(t)
        speed_list.append(v * 3.6)
        disp_list.append(x)

    return np.array(time_list), np.array(speed_list), np.array(disp_list)

# ================== 自訂 JSON 渲染（淺藍色高光）==================
LIGHT_BLUE = "#87CEEB"  # 淺天藍

def render_json_with_diff(data, default_data):
    def _format_value(value, default_value):
        is_changed = (value != default_value)
        if isinstance(value, str):
            color = LIGHT_BLUE if is_changed else "green"
            return f'<span style="color:{color};">"{value}"</span>'
        elif isinstance(value, (int, float)):
            color = LIGHT_BLUE if is_changed else "orange"
            return f'<span style="color:{color};">{value}</span>'
        else:
            return str(value)

    lines = []
    lines.append("{")
    keys = list(data.keys())
    for i, key in enumerate(keys):
        value = data[key]
        default_value = default_data.get(key)
        lines.append(f'  <span style="color:white;">"{key}"</span>: {_format_value(value, default_value)}' + ("," if i < len(keys)-1 else ""))
    lines.append("}")
    return "<br>".join(lines)

def render_battery_with_diff(battery_spec, default_battery_spec):
    if default_battery_spec is None:
        default_battery_spec = battery_spec
    lines = []
    descriptions = {
        '類型': '常見的電動載具電池類型，此處為鋰離子電池。',
        '標稱電壓 (V)': '電池組的額定電壓，由串聯電池芯數決定（每芯 3.7V）。',
        '容量 (Ah)': '電池組的總電荷容量，並聯電池芯數 × 單芯容量 (2.5Ah)。',
        '能量 (kWh)': '電池組儲存的總電能 = 電壓 × 容量 / 1000。',
        '放電倍率 (C)': '表示電池持續放電電流相對於容量的倍率，1C 代表可持續 1 小時放完電。',
        '串聯數': f'將多顆電池芯串聯以提高電壓。例如 {battery_spec["串聯數"]} 串 × 3.7V ≈ {battery_spec["串聯數"]*3.7:.0f}V。',
        '並聯數': '將多組串聯電池並聯以提高容量。總容量 = 並聯數 × 單芯容量 (2.5Ah)。',
        '估計重量 (kg)': '基於能量密度 150 Wh/kg 估算的電池組重量。'
    }
    for key in battery_spec.keys():
        value = battery_spec[key]
        default_value = default_battery_spec.get(key, value)
        is_changed = (value != default_value)
        if isinstance(value, (int, float)):
            color = LIGHT_BLUE if is_changed else "orange"
            value_str = f'<span style="color:{color};">{value}</span>'
        elif isinstance(value, str):
            color = LIGHT_BLUE if is_changed else "green"
            value_str = f'<span style="color:{color};">"{value}"</span>'
        else:
            value_str = str(value)
        lines.append(f"- **{key}**：{value_str}")
        if key in descriptions:
            lines.append(f"  > {descriptions[key]}")
    return "\n".join(lines)

# ================== Streamlit 介面 ==================
st.set_page_config(layout="centered", page_title="電動載具動力估算 (里程估計)")

st.title("⚡ 電動載具動力系統估算 (含里程估計)")

# ---------- 側邊欄（輸入參數）----------
with st.sidebar:
    st.header("🚗 整車參數規格")
    vehicle_type = st.selectbox("車種", ['小型電動車', '電動機車', '電動三輪車', '高爾夫球車'], index=1)
    weight = st.number_input("車重 (kg, 不含電池)", min_value=50, value=98, step=10)
    load = st.number_input("載重 (kg)", min_value=0, value=63, step=10)
    total_mass = weight + load
    st.caption(f"總質量: {total_mass} kg")

    speed_kmh = st.number_input("目標最高車速 (km/h)", min_value=10, value=75, step=5)
    speed_ms = speed_kmh / 3.6

    area = st.number_input("迎風面積 (m²)", min_value=0.3, value=0.61, step=0.05, format="%.2f")

    st.subheader("輪胎規格")
    tire_width = st.number_input("胎寬 (mm)", min_value=50, value=110, step=5)
    tire_aspect = st.number_input("扁平比 (%)", min_value=30, value=70, step=5)
    rim_dia_inch = st.number_input("輪胎半徑(英吋)", min_value=8, value=12, step=1,
                                   help="此處輸入的是輪輞直徑（英寸），用於計算輪胎半徑。")
    st.caption("註：此處輸入的是輪輞直徑（英寸），即輪胎內側直徑，非輪胎外徑半徑。")

    sidewall_height_mm = tire_width * tire_aspect / 100
    rim_radius_mm = (rim_dia_inch * 25.4) / 2
    tire_radius_m = (rim_radius_mm + sidewall_height_mm) / 1000
    st.caption(f"計算輪胎半徑: {tire_radius_m:.4f} m")
    wheel_radius_m = tire_radius_m

    # ---------- 動態性能表現規格 ----------
    st.header("⚡ 動態性能表現規格")
    accel_time_full = st.number_input("0→最高車速加速時間 (秒)", min_value=1.0, value=10.0, step=0.5)
    avg_accel_full = speed_ms / accel_time_full

    accel_time_0to50 = st.number_input("0→50 km/h 加速時間 (秒)", min_value=1.0, value=5.0, step=0.5)
    speed_50_ms = 50 / 3.6
    avg_accel_50 = speed_50_ms / accel_time_0to50

    # ---------- 爬坡規格 ----------
    st.header("⛰️ 爬坡規格")
    grade_percent = st.number_input("爬坡度 (%)", min_value=0.0, value=30.0, step=0.5)
    if grade_percent > 0:
        grade_angle = math.degrees(math.atan(grade_percent / 100))
        st.caption(f"換算角度: {grade_angle:.2f}°")
    else:
        st.caption("換算角度: 0°")

    # ---------- 續航里程規格 ----------
    st.header("🔋 續航里程規格")
    use_range = st.checkbox("指定續航里程 (用於電池估算)")
    if use_range:
        desired_range = st.number_input("期望續航里程 (km)", min_value=1, value=50, step=5)
    else:
        desired_range = None

    # ---------- 輸入動力鍊規格 ----------
    st.header("🔧 輸入動力鍊 (Powertrain) 規格")

    # ----- 馬達 (Motor) 區塊 -----
    with st.expander("🔹 馬達 (Motor)", expanded=True):
        voltage_option = st.radio("系統電壓", ['自動選擇', '48V', '96V'])
        if voltage_option == '自動選擇':
            voltage = None
        else:
            voltage = int(voltage_option.replace('V', ''))

        est_mode = st.radio("估算模式", ['自動估算', '手動輸入'], index=0,
                            help="自動估算：根據目標車速計算所需功率。手動輸入：您可分別設定最大功率、最大扭矩與最高轉速。")

        if est_mode == '自動估算':
            cd = get_cd_by_vehicle(vehicle_type)
            fr = FR
            required_power, _ = calculate_power_requirements(total_mass, speed_ms, area, cd, fr)
            max_power_kw = required_power * 2
            st.info(f"⚡ 所需功率 = {required_power:.2f} kW → 最大功率 = {max_power_kw:.2f} kW")
            manual_max_power = max_power_kw
            manual_peak_torque = None
            manual_max_rpm = None
        else:
            manual_max_power = st.number_input("最大功率 (kW)", min_value=0.1, value=4.4, step=0.1)
            manual_peak_torque = st.number_input("最大扭矩 (Nm)", min_value=1.0, value=18.0, step=0.1)
            manual_max_rpm = st.number_input("最高轉速 (rpm)", min_value=100, value=9000, step=100)
            base_speed_calc = (manual_max_power * 1000 * 60) / (2 * math.pi * manual_peak_torque)
            st.caption(f"對應基速 ≈ {base_speed_calc:.0f} rpm")

        motor_eff = st.number_input("馬達效率 (%)", min_value=0.0, max_value=100.0, value=90.0, step=1.0,
                                    help="固定工作點下的馬達效率，用於里程估計，也會顯示在馬達規格中。")

    # ----- 齒輪 (Gear) 區塊 -----
    with st.expander("🔹 齒輪 (Gear)", expanded=True):
        gear_option = st.radio("減速比", ['自動估算', '手動輸入'], index=1)
        if gear_option == '手動輸入':
            gear_ratio = st.number_input("請輸入減速比", min_value=1.0, value=8.7, step=0.5)
        else:
            gear_ratio = None
        gear_eff = st.number_input("齒輪箱效率 (%)", min_value=0.0, max_value=100.0, value=95.0, step=1.0,
                                   help="齒輪箱傳動效率，用於里程估計。")

    # ----- 電池 (Battery) 區塊 -----
    with st.expander("🔹 電池 (Battery)", expanded=True):
        user_battery_energy_kwh = st.number_input("電池總能量 (kWh)", min_value=0.5, value=5.0, step=0.5,
                                                  help="用於行駛里程估計的電池總能量，請根據實際電池規格輸入。")
        battery_soc = st.number_input("電池可用 SOC (%)", min_value=0.0, max_value=100.0, value=90.0, step=5.0,
                                      help="State of Charge，剩餘電量百分比，通常為避免深度放電而保留部分電量。")

    # ---------- 里程估計參數 ----------
    st.header("📊 里程估計參數")
    avg_speed_ratio = st.slider("平均車速 / 最高車速 比例", min_value=0.3, max_value=1.0, value=0.7, step=0.05,
                                help="估計里程時使用的平均車速佔最高車速的比例")
    avg_speed_kmh = speed_kmh * avg_speed_ratio
    st.caption(f"計算平均車速: {avg_speed_kmh:.1f} km/h")

    # ---------- 即時里程估計預覽（連動）----------
    cd_preview = get_cd_by_vehicle(vehicle_type)
    fr_preview = FR
    avg_speed_ms_preview = avg_speed_kmh / 3.6
    F_roll_preview = total_mass * G * fr_preview
    F_air_preview = 0.5 * RHO * cd_preview * area * avg_speed_ms_preview**2
    F_total_preview = F_roll_preview + F_air_preview
    P_wheel_preview = F_total_preview * avg_speed_ms_preview / 1000
    P_batt_preview = P_wheel_preview / (gear_eff/100) / (motor_eff/100)
    usable_energy_preview = user_battery_energy_kwh * (battery_soc / 100)
    if P_batt_preview > 0:
        hours_preview = usable_energy_preview / P_batt_preview
    else:
        hours_preview = 0
    range_preview = avg_speed_kmh * hours_preview
    st.metric("即時估計里程", f"{range_preview:.1f} km")

    st.markdown("---")
    st.caption("修改參數後，下方結果會自動更新")

# ================== 計算核心 ==================
cd = get_cd_by_vehicle(vehicle_type)
fr = FR

if voltage is None:
    if est_mode == '自動估算':
        power_val = manual_max_power
    else:
        power_val = manual_max_power
    voltage = 48 if power_val < 20 else 96

if gear_ratio is None:
    gear_ratio = estimate_gearbox(speed_ms, wheel_radius_m)

# 計算馬達最高轉速
if est_mode == '自動估算':
    required_max_rpm = speed_ms * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio
    n_max_motor = max(required_max_rpm * 1.1, 6000)
else:
    n_max_motor = manual_max_rpm

if est_mode == '自動估算':
    motor_spec, base_speed, T_peak = estimate_motor_from_power(manual_max_power, voltage, n_max_motor, motor_eff, base_speed=3000)
    max_power_kw_used = manual_max_power
else:
    motor_spec, base_speed, T_peak = estimate_motor_from_params(manual_max_power, manual_peak_torque, voltage, n_max_motor, motor_eff)
    max_power_kw_used = manual_max_power

rated_power = max_power_kw_used / 2

# 電池估算規格（僅用於顯示，不影響里程估計的電池能量）
if desired_range:
    # 注意：此處僅為了顯示電池規格，不影響里程估計（里程估計使用 user_battery_energy_kwh）
    avg_speed_ms_display = speed_ms * 0.7
    avg_speed_kmh_display = avg_speed_ms_display * 3.6
    time_h = desired_range / avg_speed_kmh_display
    avg_power_kw = rated_power * 0.7
    battery_spec = estimate_battery(avg_power_kw, voltage, duration_h=time_h)
else:
    battery_spec = estimate_battery(rated_power, voltage, duration_h=1.0)

controller_spec = estimate_controller(max_power_kw_used, voltage)

gearbox_spec = {
    '類型': '固定減速比齒輪箱',
    '減速比': round(gear_ratio, 2),
    '效率 (%)': 95
}

# ---------- 起步扭矩需求 ----------
F_roll_start = total_mass * G * fr
F_accel_full = total_mass * avg_accel_full
F_total_start_full = F_roll_start + F_accel_full
T_wheel_start_full = F_total_start_full * wheel_radius_m
T_motor_start_full = T_wheel_start_full / (gear_ratio * ETA_DRIVE)

F_accel_50 = total_mass * avg_accel_50
F_total_start_50 = F_roll_start + F_accel_50
T_wheel_start_50 = F_total_start_50 * wheel_radius_m
T_motor_start_50 = T_wheel_start_50 / (gear_ratio * ETA_DRIVE)

# ---------- 負載線 ----------
motor_rpm_flat, torque_flat,车速_flat, force_flat = calculate_load_curve(
    total_mass, area, cd, fr, wheel_radius_m, gear_ratio, speed_ms, grade_percent=0
)

if grade_percent > 0:
    motor_rpm_climb, torque_climb,车速_climb, force_climb = calculate_load_curve(
        total_mass, area, cd, fr, wheel_radius_m, gear_ratio, speed_ms, grade_percent
    )
else:
    motor_rpm_climb, torque_climb = None, None

# ---------- 加速模擬 ----------
time_acc, speed_acc, disp_acc = simulate_acceleration(
    total_mass, area, cd, fr, wheel_radius_m, gear_ratio,
    motor_spec, base_speed, T_peak, speed_ms, dt=0.1
)

# 計算實際 0-50 加速時間
if np.any(speed_acc >= 50):
    actual_0to50 = time_acc[np.argmax(speed_acc >= 50)]
else:
    actual_0to50 = np.inf

# 計算實際 0→最高車速 加速時間
if np.any(speed_acc >= speed_kmh * 0.99):
    actual_full_time = time_acc[np.argmax(speed_acc >= speed_kmh * 0.99)]
else:
    actual_full_time = np.inf

# ========== 行駛里程估計（使用側邊欄的 avg_speed_kmh，且不覆蓋）==========
# 注意：此處 avg_speed_kmh 直接來自側邊欄的滑桿計算，不受 desired_range 影響
avg_speed_ms_est = avg_speed_kmh / 3.6
F_roll_avg = total_mass * G * fr
F_air_avg = 0.5 * RHO * cd * area * avg_speed_ms_est**2
F_total_avg = F_roll_avg + F_air_avg
P_wheel_avg = F_total_avg * avg_speed_ms_est / 1000
P_battery_avg = P_wheel_avg / (gear_eff/100) / (motor_eff/100)
usable_energy = user_battery_energy_kwh * (battery_soc / 100)
if P_battery_avg > 0:
    driving_hours = usable_energy / P_battery_avg
else:
    driving_hours = 0
estimated_range = avg_speed_kmh * driving_hours

# ================== 儲存初始規格（用於比較變化）==================
if "default_motor_spec" not in st.session_state:
    st.session_state.default_motor_spec = motor_spec.copy()
if "default_battery_spec" not in st.session_state:
    st.session_state.default_battery_spec = battery_spec.copy()

# ================== 顯示區 (單欄垂直排列) ==================

# ---------- 規格摘要 ----------
st.subheader("📋 規格摘要")

with st.expander("📦 馬達規格", expanded=True):
    html = render_json_with_diff(motor_spec, st.session_state.default_motor_spec)
    st.markdown(html, unsafe_allow_html=True)

# 動態性能表現輸出
with st.expander("🏎️ 動態性能表現", expanded=True):
    st.markdown("**起步性能比較**")
    st.metric("0→最高車速起步所需馬達扭矩", f"{T_motor_start_full:.1f} Nm")
    st.metric("0→50 km/h 起步所需馬達扭矩", f"{T_motor_start_50:.1f} Nm")
    if T_motor_start_full <= T_peak and T_motor_start_50 <= T_peak:
        st.success("✅ 馬達峰值扭矩足夠滿足兩種加速需求")
    else:
        short = max(0, T_motor_start_full - T_peak, T_motor_start_50 - T_peak)
        st.error(f"❌ 馬達峰值扭矩不足，需增加 {short:.1f} Nm")
    
    st.markdown("""
    **📌 起步所需扭矩的意義：**  
    
    - 車輛從靜止開始加速時，馬達必須輸出的最小扭矩（已考慮減速比與傳動效率）。  
    - 計算時只考慮 **滾動阻力**（輪胎與地面摩擦）和 **加速慣性力** \( F = m \cdot a \)，忽略起步時的風阻。  
    - 數值越低，起步越輕快；若超過馬達峰值扭矩，則起步無力。  
    
    **❓ 為何兩個加速目標的起步扭矩不同？**  
    
    - 因為加速時間不同，導致所需的 **平均加速度** 不同。  
    - 例如：0→50 km/h 僅需 5 秒 → 平均加速度較大 → 需要的起步扭矩較高。  
    - 0→最高車速 需 10 秒 → 平均加速度較小 → 起步扭矩較低。  
    - 加速時間越短，起步瞬間需要越大的力量才能達成目標。  
    """)
    
    st.markdown("---")
    st.markdown("**加速性能對比**")
    col_a, col_b = st.columns(2)
    with col_a:
        st.metric("目標 0→50 km/h", f"{accel_time_0to50:.1f} s")
        actual_color = "green" if actual_0to50 <= accel_time_0to50 else "red"
        st.markdown(f'<p style="color:{actual_color};">實際 0→50 km/h: {actual_0to50:.1f} s</p>', unsafe_allow_html=True)
        if actual_0to50 <= accel_time_0to50:
            st.success("✅ 滿足目標")
        else:
            st.error("❌ 未達目標")
    with col_b:
        st.metric("目標 0→最高車速", f"{accel_time_full:.1f} s")
        actual_color_full = "green" if actual_full_time <= accel_time_full else "red"
        st.markdown(f'<p style="color:{actual_color_full};">實際 0→最高車速: {actual_full_time:.1f} s</p>', unsafe_allow_html=True)
        if actual_full_time <= accel_time_full:
            st.success("✅ 滿足目標")
        else:
            st.error("❌ 未達目標")
    st.caption(
        "📌 **加速性能對比說明**：\n\n"
        "• 目標值：您期望的加速時間。\n"
        "• 實際值：根據馬達最大扭矩曲線與車輛阻力模擬出的真實加速時間。\n"
        "• 若實際值 ≤ 目標值，代表馬達性能足夠；反之則需提高馬達功率或降低車重。"
    )

# 電池估算規格
with st.expander("🔋 電池 (估算規格)", expanded=False):
    battery_html = render_battery_with_diff(battery_spec, st.session_state.default_battery_spec)
    st.markdown(battery_html)

with st.expander("🎛️ 控制器", expanded=False):
    st.json(controller_spec)

with st.expander("⚙️ 齒輪箱", expanded=False):
    st.json(gearbox_spec)

with st.expander("🔁 轉換係數", expanded=False):
    torque_factor = gear_ratio * ETA_DRIVE
    speed_factor = (2 * math.pi * wheel_radius_m * 60) / (gear_ratio * 1000) * 3.6
    st.metric("輪上扭矩 / 馬達扭矩", f"{torque_factor:.3f}")
    st.caption("計算式：減速比 × 傳動效率")
    st.metric("車速 (km/h) / 馬達轉速 (rpm)", f"{speed_factor:.6f}")
    st.caption("計算式：(2π × 輪胎半徑(m) × 60) / (減速比 × 1000) × 3.6")

idx_design_local = np.argmin(np.abs(车速_flat - speed_kmh))
T_design_flat_local = (force_flat[idx_design_local] * wheel_radius_m)
F_design_flat_local = force_flat[idx_design_local]
with st.expander("🔧 設計最高車速點性能", expanded=False):
    st.metric("最高車速點輪上扭矩", f"{T_design_flat_local:.1f} Nm")
    st.metric("最高車速點輪上推力", f"{F_design_flat_local:.1f} N")

# ---------- 行駛里程估計結果展示 ----------
with st.expander("🔋 行駛里程估計", expanded=True):
    st.markdown("**估計結果**")
    col1, col2 = st.columns(2)
    with col1:
        st.metric("估計行駛里程", f"{estimated_range:.1f} km")
        st.metric("可行駛時間", f"{driving_hours:.1f} h")
    if use_range and desired_range is not None:
        with col2:
            st.metric("期望續航里程", f"{desired_range:.1f} km")
            range_met = estimated_range >= desired_range
            range_color = "green" if range_met else "red"
            st.markdown(f'<p style="color:{range_color};">比較結果: 估計里程 {"≥" if range_met else "<"} 期望里程</p>', unsafe_allow_html=True)
            if range_met:
                st.success("✅ 滿足目標")
            else:
                st.error("❌ 未達目標")
    st.markdown("---")
    st.markdown("**計算公式**")
    st.latex(r"P_{\text{avg}} = \frac{F_{\text{roll}} + F_{\text{air}}}{1000} \cdot v_{\text{avg}}")
    st.latex(r"P_{\text{batt}} = \frac{P_{\text{avg}}}{\eta_{\text{gear}} \cdot \eta_{\text{motor}}}")
    st.latex(r"E_{\text{usable}} = E_{\text{batt}} \cdot \text{SOC}")
    st.latex(r"t_{\text{drive}} = \frac{E_{\text{usable}}}{P_{\text{batt}}}")
    st.latex(r"\text{Range} = v_{\text{avg}} \cdot t_{\text{drive}}")
    st.markdown(rf"""
    - **平均車速** \(v_{{\text{{avg}}}}\) = {avg_speed_kmh:.1f} km/h  
    - **滾動阻力** \(F_{{\text{{roll}}}}\) = {F_roll_avg:.1f} N  
    - **空氣阻力** \(F_{{\text{{air}}}}\) = {F_air_avg:.1f} N  
    - **平均輪上功率** \(P_{{\text{{avg}}}}\) = {P_wheel_avg:.2f} kW  
    - **齒輪效率** \(\eta_{{\text{{gear}}}}\) = {gear_eff}%  
    - **馬達效率** \(\eta_{{\text{{motor}}}}\) = {motor_eff}%  
    - **電池平均輸出功率** \(P_{{\text{{batt}}}}\) = {P_battery_avg:.2f} kW  
    - **電池總能量 (輸入)** \(E_{{\text{{batt}}}}\) = {user_battery_energy_kwh:.2f} kWh  
    - **可用 SOC** = {battery_soc}% → **可用能量** \(E_{{\text{{usable}}}}\) = {usable_energy:.2f} kWh  
    - **可行駛時間** \(t_{{\text{{drive}}}}\) = {driving_hours:.1f} h  
    - **估計里程** = {estimated_range:.1f} km
    """)

# 下載 Excel
df_motor = pd.DataFrame([motor_spec])
df_battery = pd.DataFrame([battery_spec])
df_controller = pd.DataFrame([controller_spec])
df_gearbox = pd.DataFrame([gearbox_spec])

output = BytesIO()
with pd.ExcelWriter(output, engine='openpyxl') as writer:
    df_motor.to_excel(writer, sheet_name='馬達', index=False)
    df_battery.to_excel(writer, sheet_name='電池', index=False)
    df_controller.to_excel(writer, sheet_name='控制器', index=False)
    df_gearbox.to_excel(writer, sheet_name='齒輪箱', index=False)
st.download_button(
    label="📥 下載 Excel 報表",
    data=output.getvalue(),
    file_name="powertrain_spec.xlsx",
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    use_container_width=True
)

st.markdown("---")

# ---------- 圖1：馬達 TN 曲線 + 功率曲線 ----------
n_max_motor = motor_spec['最高轉速 (rpm)']
P_peak = motor_spec['最大功率 (kW)']

max_rpm_load = motor_rpm_flat.max() if len(motor_rpm_flat) > 0 else 0
if motor_rpm_climb is not None:
    max_rpm_load = max(max_rpm_load, motor_rpm_climb.max())
x_upper = max(n_max_motor, max_rpm_load)

n = np.linspace(0, x_upper, 500)
T_motor_max = np.zeros_like(n)
P_motor_out = np.zeros_like(n)

const_idx = n <= base_speed
T_motor_max[const_idx] = T_peak
P_motor_out[const_idx] = T_peak * n[const_idx] / 9550

power_idx = (n > base_speed) & (n <= n_max_motor)
T_motor_max[power_idx] = (P_peak * 1000) / (2 * math.pi * n[power_idx] / 60)
P_motor_out[power_idx] = P_peak

over_idx = n > n_max_motor
T_motor_max[over_idx] = 0
P_motor_out[over_idx] = 0

design_rpm = speed_ms * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio

fig1 = make_subplots(specs=[[{"secondary_y": True}]])
fig1.add_trace(
    go.Scatter(x=n, y=T_motor_max, mode='lines', name='馬達最大扭矩', line=dict(color='blue', width=3)),
    secondary_y=False
)
fig1.add_trace(
    go.Scatter(x=motor_rpm_flat, y=torque_flat, mode='lines', name='平路負載線 (馬達側扭矩)', line=dict(color='red', width=3, dash='dash')),
    secondary_y=False
)
if motor_rpm_climb is not None:
    fig1.add_trace(
        go.Scatter(x=motor_rpm_climb, y=torque_climb, mode='lines', name=f'爬坡負載線 ({grade_percent}%)',
                   line=dict(color='green', width=3, dash='dot')),
        secondary_y=False
    )
fig1.add_trace(
    go.Scatter(x=n, y=P_motor_out, mode='lines', name='馬達功率', line=dict(color='gold', width=2, dash='solid')),
    secondary_y=True
)

# 標註關鍵點
fig1.add_trace(
    go.Scatter(x=[0], y=[T_peak], mode='markers+text', name='最大扭矩點',
               text=[f'{T_peak:.1f} Nm'], textposition='bottom right',
               marker=dict(color='blue', size=12), textfont=dict(size=11)),
    secondary_y=False
)
fig1.add_trace(
    go.Scatter(x=[base_speed], y=[T_peak], mode='markers+text', name='基速點',
               text=[f'基速: {base_speed:.0f} rpm'], textposition='top left',
               marker=dict(color='green', size=12), textfont=dict(size=11)),
    secondary_y=False
)
T_at_max_n = (P_peak * 1000) / (2 * math.pi * n_max_motor / 60) if n_max_motor > 0 else 0
fig1.add_trace(
    go.Scatter(x=[n_max_motor], y=[T_at_max_n], mode='markers+text', name='最高轉速點',
               text=[f'{n_max_motor:.0f} rpm, {T_at_max_n:.1f} Nm'],
               textposition='top right',
               marker=dict(color='purple', size=12), textfont=dict(size=11)),
    secondary_y=False
)
fig1.add_vline(x=design_rpm, line_width=2, line_dash="dash", line_color="orange", opacity=0.9)
T_at_design = np.interp(design_rpm, n, T_motor_max) if design_rpm <= n_max_motor else 0
fig1.add_trace(
    go.Scatter(x=[design_rpm], y=[T_at_design], mode='markers+text',
               name='設計車速對應轉速',
               text=[f'{design_rpm:.0f} rpm, {T_at_design:.1f} Nm'],
               textposition='top center',
               marker=dict(color='orange', size=12),
               textfont=dict(size=11)),
    secondary_y=False
)

# 交點標記
intersections_flat = find_intersection(n, T_motor_max, motor_rpm_flat, torque_flat)
for i, (x_cross, y_cross) in enumerate(intersections_flat):
    fig1.add_trace(
        go.Scatter(x=[x_cross], y=[y_cross], mode='markers',
                   name=f'平路交點{i+1}' if i==0 else None,
                   marker=dict(color='red', size=12, symbol='x'),
                   showlegend=(i==0)),
        secondary_y=False
    )
    fig1.add_annotation(x=x_cross, y=y_cross,
                        text=f'{x_cross:.0f} rpm, {y_cross:.1f} Nm',
                        showarrow=True, arrowhead=2, ax=30, ay=-40,
                        font=dict(size=10))

if motor_rpm_climb is not None:
    intersections_climb = find_intersection(n, T_motor_max, motor_rpm_climb, torque_climb)
    for i, (x_cross, y_cross) in enumerate(intersections_climb):
        fig1.add_trace(
            go.Scatter(x=[x_cross], y=[y_cross], mode='markers',
                       name=f'爬坡交點{i+1}' if i==0 else None,
                       marker=dict(color='green', size=12, symbol='x'),
                       showlegend=(i==0)),
            secondary_y=False
        )
        fig1.add_annotation(x=x_cross, y=y_cross,
                            text=f'{x_cross:.0f} rpm, {y_cross:.1f} Nm',
                            showarrow=True, arrowhead=2, ax=30, ay=40,
                            font=dict(size=10))

fig1.update_layout(
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
    margin=dict(l=20, r=20, t=40, b=20),
    height=400
)
fig1.update_xaxes(title_text="轉速 (rpm)")
fig1.update_yaxes(title_text="扭矩 (Nm)", secondary_y=False)
fig1.update_yaxes(title_text="功率 (kW)", secondary_y=True)
st.plotly_chart(fig1, use_container_width=True)

st.markdown("---")

# ---------- 圖2：車輪扭矩 vs 車速 ----------
v_from_n = n / gear_ratio * (2 * math.pi * wheel_radius_m) * 3.6 / 60
T_wheel_max = T_motor_max * gear_ratio * ETA_DRIVE
T_wheel_flat = force_flat * wheel_radius_m

if grade_percent > 0:
    T_wheel_climb = force_climb * wheel_radius_m
    v_climb =车速_climb
else:
    T_wheel_climb = None
    v_climb = None

idx_design = np.argmin(np.abs(车速_flat - speed_kmh))
T_design_flat = T_wheel_flat[idx_design]
v_max_motor = n_max_motor / gear_ratio * (2 * math.pi * wheel_radius_m) * 3.6 / 60
T_at_vmax = np.interp(v_max_motor, v_from_n, T_wheel_max) if v_max_motor <= v_from_n.max() else 0

fig2 = go.Figure()
fig2.add_trace(go.Scatter(x=v_from_n, y=T_wheel_max, mode='lines', name='最大車輪扭矩',
                           line=dict(color='blue', width=3)))
fig2.add_trace(go.Scatter(x=车速_flat, y=T_wheel_flat, mode='lines', name='平路負載線',
                           line=dict(color='red', width=3, dash='dash')))
if T_wheel_climb is not None:
    fig2.add_trace(go.Scatter(x=v_climb, y=T_wheel_climb, mode='lines',
                               name=f'爬坡負載線 ({grade_percent}%)',
                               line=dict(color='green', width=3, dash='dot')))

fig2.add_vline(x=speed_kmh, line_width=2, line_dash="dash", line_color="orange", opacity=0.9)
fig2.add_trace(go.Scatter(x=[speed_kmh], y=[T_design_flat], mode='markers+text',
                           name='設計最高車速點',
                           text=[f'{speed_kmh:.0f} km/h, {T_design_flat:.1f} Nm'],
                           textposition='top right',
                           marker=dict(color='orange', size=12),
                           textfont=dict(size=11)))
fig2.add_trace(go.Scatter(x=[v_max_motor], y=[T_at_vmax], mode='markers+text',
                           name='馬達最高轉速對應車速',
                           text=[f'馬達最高速\n{v_max_motor:.0f} km/h, {T_at_vmax:.1f} Nm'],
                           textposition='top left',
                           marker=dict(color='purple', size=12),
                           textfont=dict(size=11)))

# 平路交點
intersections_flat_wheel = find_intersection(v_from_n, T_wheel_max, 车速_flat, T_wheel_flat)
for i, (x_cross, y_cross) in enumerate(intersections_flat_wheel):
    fig2.add_trace(
        go.Scatter(x=[x_cross], y=[y_cross], mode='markers',
                   name=f'平路交點{i+1}' if i==0 else None,
                   marker=dict(color='red', size=14, symbol='x'),
                   showlegend=(i==0))
    )
    fig2.add_annotation(x=x_cross, y=y_cross,
                        text=f'{x_cross:.1f} km/h',
                        showarrow=True, arrowhead=2, ax=30, ay=-40,
                        font=dict(size=10))

if T_wheel_climb is not None:
    intersections_climb_wheel = find_intersection(v_from_n, T_wheel_max, v_climb, T_wheel_climb)
    for i, (x_cross, y_cross) in enumerate(intersections_climb_wheel):
        fig2.add_trace(
            go.Scatter(x=[x_cross], y=[y_cross], mode='markers',
                       name=f'爬坡交點{i+1}' if i==0 else None,
                       marker=dict(color='green', size=14, symbol='x'),
                       showlegend=(i==0))
        )
        fig2.add_annotation(x=x_cross, y=y_cross,
                            text=f'{x_cross:.1f} km/h',
                            showarrow=True, arrowhead=2, ax=30, ay=40,
                            font=dict(size=10))
else:
    intersections_climb_wheel = []

x_vals = [speed_kmh * 1.2, v_max_motor * 1.1]
if intersections_flat_wheel:
    x_vals.extend([p[0] for p in intersections_flat_wheel])
if intersections_climb_wheel:
    x_vals.extend([p[0] for p in intersections_climb_wheel])
x_max = max(x_vals)

fig2.update_layout(
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
    margin=dict(l=20, r=20, t=40, b=20),
    height=400
)
fig2.update_xaxes(title_text="車速 (km/h)", range=[0, x_max])
fig2.update_yaxes(title_text="扭矩 (Nm)")
st.plotly_chart(fig2, use_container_width=True)

st.markdown("---")

# ---------- 圖3：速度與位移 vs 時間 ----------
fig3 = make_subplots(specs=[[{"secondary_y": True}]])
fig3.add_trace(
    go.Scatter(x=time_acc, y=speed_acc, mode='lines', name='車速 (km/h)', line=dict(color='blue', width=3)),
    secondary_y=False
)
fig3.add_trace(
    go.Scatter(x=time_acc, y=disp_acc, mode='lines', name='位移 (m)', line=dict(color='red', width=2, dash='dash')),
    secondary_y=True
)

# 標註達到 50 km/h 的時間
idx_50 = np.argmax(speed_acc >= 50)
if idx_50 > 0:
    t_50 = time_acc[idx_50]
    fig3.add_vline(x=t_50, line_width=1, line_dash="dot", line_color="orange", opacity=0.7)
    fig3.add_annotation(x=t_50, y=50, text=f"實際50km/h @ {t_50:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=-30)

# 標註達到最高車速的時間
idx_max = np.argmax(speed_acc >= speed_kmh * 0.99)
if idx_max > 0:
    t_max = time_acc[idx_max]
    fig3.add_vline(x=t_max, line_width=1, line_dash="dot", line_color="green", opacity=0.7)
    fig3.add_annotation(x=t_max, y=speed_kmh, text=f"實際{int(speed_kmh)}km/h @ {t_max:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=30)

# 目標 0-50 km/h 加速時間紫色虛線
fig3.add_vline(x=accel_time_0to50, line_width=1, line_dash="dot", line_color="purple", opacity=0.7)
fig3.add_annotation(x=accel_time_0to50, y=50, text=f"目標50km/h @ {accel_time_0to50:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=-50, font=dict(color="purple"))

# 目標 0→最高車速 加速時間棕色虛線
fig3.add_vline(x=accel_time_full, line_width=1, line_dash="dot", line_color="brown", opacity=0.7)
fig3.add_annotation(x=accel_time_full, y=speed_kmh, text=f"目標最高車速 @ {accel_time_full:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=-70, font=dict(color="brown"))

fig3.update_layout(
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
    margin=dict(l=20, r=20, t=40, b=20),
    height=400
)
fig3.update_xaxes(title_text="時間 (秒)")
fig3.update_yaxes(title_text="車速 (km/h)", secondary_y=False)
fig3.update_yaxes(title_text="位移 (m)", secondary_y=True)
st.plotly_chart(fig3, use_container_width=True)

st.markdown("---")
st.caption("💡 提示：圖中紫色虛線為目標 0→50 km/h 加速時間，棕色虛線為目標 0→最高車速加速時間。")