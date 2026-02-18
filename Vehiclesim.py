import streamlit as st
import pandas as pd
import numpy as np
import math
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from io import BytesIO

# ================== 常數與預設參數 ==================
G = 9.81
RHO = 1.2
FR = 0.015
CD = 0.4
ETA_DRIVE = 0.9                # 傳動效率
ETA_MOTOR = 1.0                 # 馬達效率（設為1，不影響計算）
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

def estimate_motor_from_power(max_power_kw, voltage, n_max, base_speed=3000):
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
        '效率 (%)': round(ETA_MOTOR * 100, 1),
        '估計重量 (kg)': round(max_power_kw * MOTOR_POWER_DENSITY, 1)
    }
    return motor_spec, base_speed, T_peak

def estimate_motor_from_params(max_power_kw, peak_torque_Nm, voltage, n_max):
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
        '效率 (%)': round(ETA_MOTOR * 100, 1),
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
        F_wheel = T_motor * gear_ratio * ETA_DRIVE / wheel_radius_m
        F_resist = resistance(v)
        F_net = F_wheel - F_resist
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

# ================== Streamlit 介面 ==================
st.set_page_config(layout="centered", page_title="電動載具動力估算 v1.0 優化版")

st.title("⚡ 電動載具動力系統估算 (優化版 v1.0)")

# ---------- 側邊欄（輸入參數）----------
with st.sidebar:
    st.header("🚗 輸入參數")
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

    # 系統電壓預設為 48V
    voltage_option = st.radio("系統電壓", ['自動選擇', '48V', '96V'], index=1)  # 預設 48V
    if voltage_option == '自動選擇':
        voltage = None
    else:
        voltage = int(voltage_option.replace('V', ''))

    # 減速比預設為手動輸入，值 8.7
    gear_option = st.radio("減速比", ['自動估算', '手動輸入'], index=1)  # 預設手動輸入
    if gear_option == '手動輸入':
        gear_ratio = st.number_input("請輸入減速比", min_value=1.0, value=8.7, step=0.5)
    else:
        gear_ratio = None

    # ---------- 馬達規格預估 ----------
    st.markdown("---")
    st.subheader("🔧 馬達規格預估")
    # 預設為手動輸入
    est_mode = st.radio("估算模式", ['自動估算', '手動輸入'], index=1,
                        help="自動估算：根據目標車速計算所需功率。手動輸入：您可分別設定最大功率與最大扭矩。")

    if est_mode == '自動估算':
        cd = get_cd_by_vehicle(vehicle_type)
        fr = FR
        required_power, _ = calculate_power_requirements(total_mass, speed_ms, area, cd, fr)
        max_power_kw = required_power * 2
        st.info(f"⚡ 所需功率 = {required_power:.2f} kW → 最大功率 = {max_power_kw:.2f} kW")
        manual_max_power = max_power_kw
        manual_peak_torque = None
    else:
        # 手動輸入預設值 4.2 kW, 17 Nm
        manual_max_power = st.number_input("最大功率 (kW)", min_value=0.1, value=4.2, step=0.1)
        manual_peak_torque = st.number_input("最大扭矩 (Nm)", min_value=1.0, value=17.0, step=0.1)
        base_speed_calc = (manual_max_power * 1000 * 60) / (2 * math.pi * manual_peak_torque)
        st.caption(f"對應基速 ≈ {base_speed_calc:.0f} rpm")

    # ---------- 加速度規格 ----------
    st.markdown("---")
    st.subheader("⚡ 加速度規格")
    # 預設 15 秒
    accel_time_full = st.number_input("0→最高車速加速時間 (秒)", min_value=1.0, value=15.0, step=0.5)
    avg_accel_full = speed_ms / accel_time_full

    accel_time_0to50 = st.number_input("0→50 km/h 加速時間 (秒)", min_value=1.0, value=5.0, step=0.5)
    speed_50_ms = 50 / 3.6
    avg_accel_50 = speed_50_ms / accel_time_0to50

    # ---------- 爬坡度設定 ----------
    st.markdown("---")
    st.subheader("⛰️ 爬坡設定")
    # 預設 18%
    grade_percent = st.number_input("爬坡度 (%)", min_value=0.0, value=18.0, step=0.5)

    # ---------- 續航里程設定 ----------
    st.markdown("---")
    st.subheader("🔋 續航設定")
    use_range = st.checkbox("指定續航里程 (用於電池估算)", value=False)
    if use_range:
        desired_range = st.number_input("期望續航里程 (km)", min_value=1, value=50, step=5)
    else:
        desired_range = None

    st.markdown("---")
    st.caption("修改參數後，下方結果會自動更新")

# ================== 計算核心 ==================
cd = get_cd_by_vehicle(vehicle_type)
fr = FR

# 確定電壓（若自動選擇則根據功率）
if voltage is None:
    if est_mode == '自動估算':
        power_val = manual_max_power
    else:
        power_val = manual_max_power
    voltage = 48 if power_val < 20 else 96

# 確定減速比
if gear_ratio is None:
    gear_ratio = estimate_gearbox(speed_ms, wheel_radius_m)

# 計算馬達最高轉速
required_max_rpm = speed_ms * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio
n_max_motor = max(required_max_rpm * 1.1, 6000)

# 根據模式估算馬達
if est_mode == '自動估算':
    motor_spec, base_speed, T_peak = estimate_motor_from_power(manual_max_power, voltage, n_max_motor, base_speed=3000)
    max_power_kw_used = manual_max_power
else:
    motor_spec, base_speed, T_peak = estimate_motor_from_params(manual_max_power, manual_peak_torque, voltage, n_max_motor)
    max_power_kw_used = manual_max_power

rated_power = max_power_kw_used / 2

# 電池估算（根據續航選項）
if desired_range:
    # 假設平均車速為最高車速的 0.7，平均功率為額定功率的 0.7
    avg_speed_ms = speed_ms * 0.7
    avg_speed_kmh = avg_speed_ms * 3.6
    time_h = desired_range / avg_speed_kmh
    avg_power_kw = rated_power * 0.7
    battery_spec = estimate_battery(avg_power_kw, voltage, duration_h=time_h)
else:
    # 預設以額定功率運行 1 小時
    battery_spec = estimate_battery(rated_power, voltage, duration_h=1.0)

# 控制器
controller_spec = estimate_controller(max_power_kw_used, voltage)

# 齒輪箱
gearbox_spec = {
    '類型': '固定減速比齒輪箱',
    '減速比': round(gear_ratio, 2),
    '效率 (%)': 95
}

# ---------- 定速 30 km/h 續航需求計算 ----------
CRUISE_SPEED_KMH = 30.0
CRUISE_RANGE_KM = 90.0
cruise_speed_ms = CRUISE_SPEED_KMH / 3.6
# 計算定速阻力
F_roll_cruise = total_mass * G * fr
F_air_cruise = 0.5 * RHO * cd * area * cruise_speed_ms**2
F_total_cruise = F_roll_cruise + F_air_cruise
P_wheel_cruise = F_total_cruise * cruise_speed_ms
P_motor_cruise = P_wheel_cruise / ETA_DRIVE / 1000  # kW
# 所需時間
time_h_cruise = CRUISE_RANGE_KM / CRUISE_SPEED_KMH
# 所需能量 (Wh)
energy_needed_wh = P_motor_cruise * 1000 * time_h_cruise
# 對應的電池容量 (Ah) at 系統電壓
capacity_needed_ah = energy_needed_wh / voltage
# 當前電池能量 (Wh)
current_battery_wh = battery_spec['能量 (kWh)'] * 1000
battery_enough = current_battery_wh >= energy_needed_wh

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

# ---------- 0-50加速需求曲線 ----------
F_accel_const = total_mass * avg_accel_50
T_accel_const_motor = F_accel_const * wheel_radius_m / (gear_ratio * ETA_DRIVE)
v_accel = np.linspace(0, min(50, speed_kmh), 50)
torque_flat_at_v = np.interp(v_accel, 车速_flat, torque_flat)
torque_total_accel_motor = torque_flat_at_v + T_accel_const_motor
torque_total_accel_wheel = torque_total_accel_motor * gear_ratio * ETA_DRIVE

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

# ================== 顯示區 (單欄垂直排列) ==================

# ---------- 規格摘要 ----------
st.subheader("📋 規格摘要")

with st.expander("📦 馬達規格", expanded=True):
    st.json(motor_spec)

with st.expander("🚦 起步性能比較", expanded=True):
    st.metric("0→最高車速起步所需馬達扭矩", f"{T_motor_start_full:.1f} Nm")
    st.metric("0→50 km/h 起步所需馬達扭矩", f"{T_motor_start_50:.1f} Nm")
    if T_motor_start_full <= T_peak and T_motor_start_50 <= T_peak:
        st.success("✅ 馬達峰值扭矩足夠滿足兩種加速需求")
    else:
        short = max(0, T_motor_start_full - T_peak, T_motor_start_50 - T_peak)
        st.error(f"❌ 馬達峰值扭矩不足，需增加 {short:.1f} Nm")

with st.expander("⚡ 加速性能對比", expanded=True):
    col_a, col_b = st.columns(2)
    with col_a:
        st.metric("目標 0→50 km/h", f"{accel_time_0to50:.1f} s")
        st.metric("實際 0→50 km/h", f"{actual_0to50:.1f} s")
        if actual_0to50 <= accel_time_0to50:
            st.success("✅ 滿足目標")
        else:
            st.error("❌ 未達目標")
    with col_b:
        st.metric("目標 0→最高車速", f"{accel_time_full:.1f} s")
        st.metric("實際 0→最高車速", f"{actual_full_time:.1f} s")
        if actual_full_time <= accel_time_full:
            st.success("✅ 滿足目標")
        else:
            st.error("❌ 未達目標")

with st.expander("🔋 電池", expanded=False):
    st.json(battery_spec)
    st.markdown("---")
    st.markdown(f"**定速 {CRUISE_SPEED_KMH} km/h 行駛 {CRUISE_RANGE_KM} km 需求**")
    st.metric("所需能量", f"{energy_needed_wh:.0f} Wh")
    st.metric("所需容量 (@{voltage}V)", f"{capacity_needed_ah:.1f} Ah")
    if battery_enough:
        st.success("✅ 當前電池能量足夠")
    else:
        st.error(f"❌ 當前電池能量不足，短缺 {energy_needed_wh - current_battery_wh:.0f} Wh")

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
               marker=dict(color='blue', size=10), textfont=dict(size=10)),
    secondary_y=False
)
fig1.add_trace(
    go.Scatter(x=[base_speed], y=[T_peak], mode='markers+text', name='基速點',
               text=[f'基速: {base_speed:.0f} rpm'], textposition='top left',
               marker=dict(color='green', size=10), textfont=dict(size=10)),
    secondary_y=False
)
T_at_max_n = (P_peak * 1000) / (2 * math.pi * n_max_motor / 60) if n_max_motor > 0 else 0
fig1.add_trace(
    go.Scatter(x=[n_max_motor], y=[T_at_max_n], mode='markers+text', name='最高轉速點',
               text=[f'{n_max_motor:.0f} rpm, {T_at_max_n:.1f} Nm'],
               textposition='top right',
               marker=dict(color='purple', size=10), textfont=dict(size=10)),
    secondary_y=False
)
fig1.add_vline(x=design_rpm, line_width=2, line_dash="dash", line_color="orange", opacity=0.9)
T_at_design = np.interp(design_rpm, n, T_motor_max) if design_rpm <= n_max_motor else 0
fig1.add_trace(
    go.Scatter(x=[design_rpm], y=[T_at_design], mode='markers+text',
               name='設計車速對應轉速',
               text=[f'{design_rpm:.0f} rpm, {T_at_design:.1f} Nm'],
               textposition='top center',
               marker=dict(color='orange', size=10),
               textfont=dict(size=10)),
    secondary_y=False
)

# 交點
intersections_flat = find_intersection(n, T_motor_max, motor_rpm_flat, torque_flat)
for i, (x_cross, y_cross) in enumerate(intersections_flat):
    fig1.add_trace(
        go.Scatter(x=[x_cross], y=[y_cross], mode='markers',
                   name=f'平路交點{i+1}' if i==0 else None,
                   marker=dict(color='red', size=10, symbol='x'),
                   showlegend=(i==0)),
        secondary_y=False
    )
    fig1.add_annotation(x=x_cross, y=y_cross,
                        text=f'{x_cross:.0f} rpm, {y_cross:.1f} Nm',
                        showarrow=True, arrowhead=2, ax=20, ay=-30,
                        font=dict(size=9))

if motor_rpm_climb is not None:
    intersections_climb = find_intersection(n, T_motor_max, motor_rpm_climb, torque_climb)
    for i, (x_cross, y_cross) in enumerate(intersections_climb):
        fig1.add_trace(
            go.Scatter(x=[x_cross], y=[y_cross], mode='markers',
                       name=f'爬坡交點{i+1}' if i==0 else None,
                       marker=dict(color='green', size=10, symbol='x'),
                       showlegend=(i==0)),
            secondary_y=False
        )
        fig1.add_annotation(x=x_cross, y=y_cross,
                            text=f'{x_cross:.0f} rpm, {y_cross:.1f} Nm',
                            showarrow=True, arrowhead=2, ax=20, ay=30,
                            font=dict(size=9))

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
fig2.add_trace(go.Scatter(x=v_accel, y=torque_total_accel_wheel, mode='lines',
                           name=f'0-50km/h加速需求 ({accel_time_0to50}s)',
                           line=dict(color='blue', width=2, dash='dash')))

fig2.add_vline(x=speed_kmh, line_width=2, line_dash="dash", line_color="orange", opacity=0.9)
fig2.add_trace(go.Scatter(x=[speed_kmh], y=[T_design_flat], mode='markers+text',
                           name='設計最高車速點',
                           text=[f'{speed_kmh:.0f} km/h, {T_design_flat:.1f} Nm'],
                           textposition='top right',
                           marker=dict(color='orange', size=10),
                           textfont=dict(size=10)))
fig2.add_trace(go.Scatter(x=[v_max_motor], y=[T_at_vmax], mode='markers+text',
                           name='馬達最高轉速對應車速',
                           text=[f'馬達最高速\n{v_max_motor:.0f} km/h, {T_at_vmax:.1f} Nm'],
                           textposition='top left',
                           marker=dict(color='purple', size=10),
                           textfont=dict(size=10)))

intersections_flat_wheel = find_intersection(v_from_n, T_wheel_max, 车速_flat, T_wheel_flat)
for i, (x_cross, y_cross) in enumerate(intersections_flat_wheel):
    fig2.add_trace(
        go.Scatter(x=[x_cross], y=[y_cross], mode='markers',
                   name=f'平路交點{i+1}' if i==0 else None,
                   marker=dict(color='red', size=12, symbol='x'),
                   showlegend=(i==0))
    )
    fig2.add_annotation(x=x_cross, y=y_cross,
                        text=f'{x_cross:.1f} km/h',
                        showarrow=True, arrowhead=2, ax=20, ay=-30,
                        font=dict(size=9))

if T_wheel_climb is not None:
    intersections_climb_wheel = find_intersection(v_from_n, T_wheel_max, v_climb, T_wheel_climb)
    for i, (x_cross, y_cross) in enumerate(intersections_climb_wheel):
        fig2.add_trace(
            go.Scatter(x=[x_cross], y=[y_cross], mode='markers',
                       name=f'爬坡交點{i+1}' if i==0 else None,
                       marker=dict(color='green', size=12, symbol='x'),
                       showlegend=(i==0))
        )
        fig2.add_annotation(x=x_cross, y=y_cross,
                            text=f'{x_cross:.1f} km/h',
                            showarrow=True, arrowhead=2, ax=20, ay=30,
                            font=dict(size=9))
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