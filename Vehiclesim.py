import streamlit as st
import pandas as pd
import numpy as np
import math
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from io import BytesIO
import re

# ================== 隱藏版本歷史記錄（不顯示於頁面）==================
"""
v3.1 (2026-04-26) - 交點數值標籤版
- 圖1與圖2中負載線與馬達曲線的交點新增數值標籤（自動偏移防重疊）。
- 統一圖1、圖2、圖5標示名稱，移除不必要的能耗傳遞鏈說明。
- 圖2改為「車輪推力」取代車輪扭矩。
- 移除圖5中的「驅動輸出能耗」顯示。

v3.0 (2026-04-25) - 效率地圖整合版
- 整合效率地圖支援（標準三欄／寬格式）。
- 新增圖5能耗積分與里程估算。
- 新增圖6效率地圖等高線與工作點，工作點超限紅標。
- 動態計算範例，靜默化讀取流程。

v2.1 (2026-04-11) - TN曲線支援版
- 新增讀取馬達TN曲線模式；加速模擬支援自訂曲線；圖表標籤防重疊優化。

v2.0 (2026-04-06) - WLTC整合版
- 整合WLTC行駛工況分析、負扭矩顯示、反向動力學計算。

v1.2 (2026-04-06) - 圖表標籤優化
- 圖表標籤優化、負載線交點重新命名、圖2 Y軸範圍調整。

v1.1 (2026-04-05) - 介面重組
- 側邊欄介面重組、負載線延伸至馬達極速、數值高光改淺藍色、新增里程比較。

v1.0 (2026-02-19) - 初始版本
- 基礎輸入、三張圖表、加速對比、Excel 下載、手機優化。
"""
# ================== 檢查 SciPy 可用性 ==================
try:
    from scipy import integrate
    from scipy.interpolate import griddata
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    st.warning("⚠️ 未安裝 SciPy，部分進階功能（效率地圖插值、累積積分）將無法使用。請執行 `pip install scipy` 以啟用完整功能。")

# ================== 常數與預設參數 ==================
G = 9.81
RHO = 1.2
FR = 0.015
CD = 0.4
ETA_DRIVE = 0.9
ETA_MOTOR = 1.0
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

def calculate_load_curve(mass, area, cd, fr, wheel_radius_m, gear_ratio, speed_max_ms, grade_percent=0, extend_to_vmax=None):
    if extend_to_vmax is not None:
        max_speed = max(speed_max_ms, extend_to_vmax)
    else:
        max_speed = speed_max_ms * 1.1
    speeds = np.linspace(0, max_speed, 150)
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

def simulate_acceleration(mass, area, cd, fr, wheel_radius_m, gear_ratio, motor_spec, base_speed, T_peak, speed_max_ms, dt=0.1, custom_tn_df=None):
    n_max = motor_spec['最高轉速 (rpm)']
    P_peak = motor_spec['最大功率 (kW)']

    def get_max_torque(v):
        if v <= 0:
            v_n = 0
        else:
            v_n = v * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio
            
        if custom_tn_df is not None:
            return np.interp(v_n, custom_tn_df['rpm'].values, custom_tn_df['torque'].values, right=0)
            
        if v_n <= base_speed:
            return T_peak
        elif v_n <= n_max:
            return (P_peak * 1000) / (2 * math.pi * v_n / 60)
        else:
            return 0

    t, v, x = 0, 0, 0
    time_list, speed_list, disp_list = [0], [0], [0]

    def resistance(v_ms):
        return (mass * G * fr) + (0.5 * RHO * cd * area * v_ms**2)

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

# ================== 行駛工況工作點計算 ==================
def compute_motor_operating_points_from_wltc(df_wltc, mass, area, cd, fr, wheel_radius_m, gear_ratio, gear_eff):
    times = df_wltc['time'].values
    speeds_kmh = df_wltc['speed_kmh'].values
    
    if 'accel_ms2' in df_wltc.columns:
        accels = df_wltc['accel_ms2'].values
    else:
        speeds_ms = speeds_kmh / 3.6
        accels = np.gradient(speeds_ms, times)
        df_wltc['accel_ms2'] = accels
    
    n = len(times)
    F_roll_const = mass * G * fr
    factor_aero = 0.5 * RHO * cd * area
    wheel_circ = 2 * math.pi * wheel_radius_m
    motor_rpm = np.zeros(n)
    motor_torque = np.zeros(n)
    wheel_torque = np.zeros(n)
    
    for i in range(n):
        v_ms = speeds_kmh[i] / 3.6
        F_roll = F_roll_const
        F_air = factor_aero * v_ms**2
        F_accel = mass * accels[i]
        F_tractive = F_roll + F_air + F_accel
        motor_rpm[i] = v_ms / wheel_circ * 60 * gear_ratio
        
        if F_tractive >= 0:
            motor_torque[i] = F_tractive * wheel_radius_m / (gear_ratio * gear_eff)
        else:
            motor_torque[i] = F_tractive * wheel_radius_m / (gear_ratio * gear_eff)
            
        wheel_torque[i] = motor_torque[i] * gear_ratio * gear_eff
        
    result_df = pd.DataFrame({
        'time': times,
        'speed_kmh': speeds_kmh,
        'accel_ms2': accels,
        'motor_rpm': motor_rpm,
        'motor_torque_Nm': motor_torque,
        'wheel_torque_Nm': wheel_torque
    })
    return result_df

# ================== 理論能耗計算函數 ==================
def build_efficiency_interpolator(df_eff):
    if df_eff is None:
        return None
    rpm = df_eff['rpm'].values
    torque = df_eff['torque'].values
    eff = df_eff['efficiency'].values
    points = np.vstack((rpm, torque)).T
    def interp_func(rpm_query, torque_query):
        if np.isscalar(rpm_query):
            rpm_query = np.array([rpm_query])
            torque_query = np.array([torque_query])
        else:
            rpm_query = np.asarray(rpm_query)
            torque_query = np.asarray(torque_query)
        query_points = np.vstack((rpm_query, torque_query)).T
        eff_interp = griddata(points, eff, query_points, method='linear', fill_value=np.nan)
        nan_mask = np.isnan(eff_interp)
        if np.any(nan_mask):
            eff_nearest = griddata(points, eff, query_points, method='nearest')
            eff_interp[nan_mask] = eff_nearest[nan_mask]
        return eff_interp
    return interp_func

def compute_theoretical_energy_consumption(df_cycle, mass, area, cd, fr, gear_eff_percent, motor_eff_percent, eff_interpolator=None, df_operating_points=None):
    times = df_cycle['time'].values
    speeds_kmh = df_cycle['speed_kmh'].values
    accels = df_cycle['accel_ms2'].values
    speeds_ms = speeds_kmh / 3.6
    
    F_roll = mass * G * fr
    F_aero = 0.5 * RHO * cd * area * speeds_ms**2
    F_acc = mass * accels
    F_trac = F_roll + F_aero + F_acc
    P_wheel = F_trac * speeds_ms
    
    if eff_interpolator is not None and df_operating_points is not None:
        motor_rpm = df_operating_points['motor_rpm'].values
        motor_torque = df_operating_points['motor_torque_Nm'].values
        eff_values = eff_interpolator(motor_rpm, motor_torque)
        eff_values = np.clip(eff_values, 1, 100) / 100.0
        eta_gear = gear_eff_percent / 100.0
        eta_total = eta_gear * eff_values
    else:
        eta_gear = gear_eff_percent / 100.0
        eta_motor_fixed = motor_eff_percent / 100.0
        eta_total = eta_gear * eta_motor_fixed
    
    P_batt = np.zeros_like(P_wheel)
    drive_mask = P_wheel >= 0
    regen_mask = P_wheel < 0
    if eff_interpolator is not None and df_operating_points is not None:
        P_batt[drive_mask] = P_wheel[drive_mask] / eta_total[drive_mask]
        P_batt[regen_mask] = P_wheel[regen_mask] * eta_total[regen_mask]
    else:
        P_batt[drive_mask] = P_wheel[drive_mask] / eta_total
        P_batt[regen_mask] = P_wheel[regen_mask] * eta_total
    
    P_drive_batt = np.maximum(P_batt, 0)
    P_regen_batt = np.minimum(P_batt, 0)
    
    if SCIPY_AVAILABLE:
        try:
            cum_trapz_func = integrate.cumulative_trapezoid
        except AttributeError:
            cum_trapz_func = integrate.cumtrapz
        energy_wh_cumulative = cum_trapz_func(P_batt, times, initial=0) / 3600.0
        total_energy_wh = energy_wh_cumulative[-1]
        drive_energy_wh = integrate.trapezoid(P_drive_batt, times) / 3600.0
        regen_energy_wh = integrate.trapezoid(P_regen_batt, times) / 3600.0
        wheel_energy_wh = integrate.trapezoid(P_wheel, times) / 3600.0
        distance_km = integrate.trapezoid(speeds_ms, times) / 1000.0
    else:
        total_energy_wh = np.trapz(P_batt, times) / 3600.0
        drive_energy_wh = np.trapz(P_drive_batt, times) / 3600.0
        regen_energy_wh = np.trapz(P_regen_batt, times) / 3600.0
        wheel_energy_wh = np.trapz(P_wheel, times) / 3600.0
        distance_km = np.trapz(speeds_ms, times) / 1000.0
        dt = np.diff(times, prepend=times[0])
        energy_wh_cumulative = np.cumsum(P_batt * dt) / 3600.0
    
    wh_per_km_batt = total_energy_wh / distance_km if distance_km > 0 else 0
    wh_per_km_wheel = wheel_energy_wh / distance_km if distance_km > 0 else 0
    
    return wh_per_km_batt, wh_per_km_wheel, distance_km, total_energy_wh, drive_energy_wh, regen_energy_wh, times, P_batt, energy_wh_cumulative

# ================== 求交點函數 ==================
def find_intersection(x1, y1, x2, y2):
    """求兩條曲線的交點 (x1,y1) 和 (x2,y2)"""
    y2_interp = np.interp(x1, x2, y2)
    diff = y1 - y2_interp
    intersections = []
    for i in range(len(x1)-1):
        if diff[i] * diff[i+1] <= 0:
            # 線性內插交點
            x_cross = x1[i] - diff[i] * (x1[i+1] - x1[i]) / (diff[i+1] - diff[i])
            y_cross = np.interp(x_cross, x1, y1)
            intersections.append((x_cross, y_cross))
    return intersections

# ================== 自訂 JSON 渲染 ==================
LIGHT_BLUE = "#87CEEB"

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

    lines = ["{"]
    keys = list(data.keys())
    for i, key in enumerate(keys):
        value = data[key]
        default_value = default_data.get(key)
        lines.append(f'  <span style="color:white;">"{key}"</span>: {_format_value(value, default_value)}' + ("," if i < len(keys)-1 else ""))
    lines.append("}")
    return "<br>".join(lines)

def render_battery_with_diff(battery_spec, default_battery_spec):
    if default_battery_spec is None: default_battery_spec = battery_spec
    lines = []
    descriptions = {
        '類型': '常見的電動載具電池類型，此處為鋰離子電池。',
        '標稱電壓 (V)': '電池組的額定電壓，由串聯電池芯數決定（每芯 3.7V）。',
        '容量 (Ah)': '電池組的總電荷容量，並聯電池芯數 × 單芯容量 (2.5Ah)。',
        '能量 (kWh)': '電池組儲存的總電能 = 電壓 × 容量 / 1000。',
        '放電倍率 (C)': '表示電池持續放電電流相對於容量的倍率，1C 代表可持续 1 小時放完電。',
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
st.set_page_config(layout="centered", page_title="電動載具動力系統估算")

st.title("⚡ 電動載具動力系統估算")

# ---------- 側邊欄（輸入參數）----------
with st.sidebar:
    st.header("🚗 整車參數規格")
    
    weight = st.number_input("車重 (kg, 不含電池)", min_value=50, value=98, step=10)
    load = st.number_input("載重 (kg)", min_value=0, value=63, step=10)
    total_mass = weight + load
    st.caption(f"總質量: {total_mass} kg")

    speed_kmh = st.number_input("目標最高車速 (km/h)", min_value=10, value=75, step=5)
    speed_ms = speed_kmh / 3.6

    st.subheader("🌬️ 阻力與環境參數設定")
    area = st.number_input("迎風面積 A (m²)", min_value=0.1, value=0.61, step=0.01, format="%.2f")
    cd = st.number_input("風阻係數 Cd", min_value=0.1, max_value=2.0, value=0.50, step=0.01, format="%.2f")
    fr = st.number_input("滾動阻力係數 fr", min_value=0.001, max_value=0.05, value=0.015, step=0.001, format="%.3f")

    st.session_state.cd = cd
    st.session_state.fr = fr

    st.subheader("輪胎規格")
    tire_width = st.number_input("胎寬 (mm)", min_value=50, value=110, step=5)
    tire_aspect = st.number_input("扁平比 (%)", min_value=30, value=70, step=5)
    rim_dia_inch = st.number_input("輪胎半徑(英吋)", min_value=8, value=12, step=1)
    sidewall_height_mm = tire_width * tire_aspect / 100
    rim_radius_mm = (rim_dia_inch * 25.4) / 2
    tire_radius_m = (rim_radius_mm + sidewall_height_mm) / 1000
    st.caption(f"計算輪胎半徑: {tire_radius_m:.4f} m")
    wheel_radius_m = tire_radius_m
    st.session_state.wheel_radius_m = wheel_radius_m

    st.header("⚡ 動態性能表現規格")
    accel_time_full = st.number_input("0→最高車速加速時間 (秒)", min_value=1.0, value=10.0, step=0.5)
    avg_accel_full = speed_ms / accel_time_full
    accel_time_0to50 = st.number_input("0→50 km/h 加速時間 (秒)", min_value=1.0, value=5.0, step=0.5)
    speed_50_ms = 50 / 3.6
    avg_accel_50 = speed_50_ms / accel_time_0to50

    st.header("⛰️ 爬坡規格")
    grade_percent = st.number_input("爬坡度 (%)", min_value=0.0, value=30.0, step=0.5)
    if grade_percent > 0:
        grade_angle = math.degrees(math.atan(grade_percent / 100))
        st.caption(f"換算角度: {grade_angle:.2f}°")
    else:
        st.caption("換算角度: 0°")

    st.header("🔋 續航里程規格")
    use_range = st.checkbox("指定續航里程 (用於電池估算)")
    if use_range:
        desired_range = st.number_input("期望續航里程 (km)", min_value=1, value=50, step=5)
    else:
        desired_range = None

    st.header("🔧 輸入動力鍊 (Powertrain) 規格")

    with st.expander("🔹 馬達 (Motor)", expanded=True):
        voltage_option = st.radio("系統電壓", ['自動選擇', '48V', '96V'])
        if voltage_option == '自動選擇':
            voltage = None
        else:
            voltage = int(voltage_option.replace('V', ''))

        est_mode = st.radio("估算模式", ['手動輸入', '自動估算', '讀取馬達TN曲線'], index=0)

        if est_mode == '自動估算':
            required_power, _ = calculate_power_requirements(total_mass, speed_ms, area, cd, fr)
            max_power_kw = required_power * 2
            st.info(f"⚡ 所需功率 = {required_power:.2f} kW → 最大功率 = {max_power_kw:.2f} kW")
            manual_max_power = max_power_kw
            manual_peak_torque = None
            manual_max_rpm = None
            custom_tn_df = None
        elif est_mode == '手動輸入':
            manual_max_power = st.number_input("最大功率 (kW)", min_value=0.1, value=4.4, step=0.1)
            manual_peak_torque = st.number_input("最大扭矩 (Nm)", min_value=1.0, value=18.0, step=0.1)
            manual_max_rpm = st.number_input("最高轉速 (rpm)", min_value=100, value=9000, step=100)
            base_speed_calc = (manual_max_power * 1000 * 60) / (2 * math.pi * manual_peak_torque)
            st.caption(f"對應基速 ≈ {base_speed_calc:.0f} rpm")
            custom_tn_df = None
        else:
            st.markdown("請上傳包含轉速(rpm)與扭矩(Nm)的 CSV 檔案")
            tn_file = st.file_uploader("上傳 TN 曲線 (CSV)", type=["csv"], key="tn_upload")
            custom_tn_df = None
            manual_max_power = 4.4 
            manual_peak_torque = 18.0
            manual_max_rpm = 9000
            if tn_file is not None:
                try:
                    df_tn = pd.read_csv(tn_file)
                    st.success("成功讀取多重 TN 曲線檔案")
                    rpm_col = st.selectbox("👉 選擇轉速 (X軸) 欄位", df_tn.columns, index=0)
                    available_torque_curves = [col for col in df_tn.columns if col != rpm_col]
                    if not available_torque_curves:
                        st.error("檔案中除了轉速外，找不到其他的扭力數據欄位！")
                    else:
                        t_col = st.selectbox("👉 選擇本次模擬使用的扭力曲線", available_torque_curves, index=0)
                        df_tn_active = df_tn.sort_values(by=rpm_col).dropna(subset=[rpm_col, t_col])
                        custom_tn_df = df_tn_active[[rpm_col, t_col]].copy()
                        custom_tn_df.columns = ['rpm', 'torque']
                        custom_tn_df['power_kw'] = custom_tn_df['torque'] * custom_tn_df['rpm'] / 9550.0
                        manual_max_rpm = custom_tn_df['rpm'].max()
                        manual_peak_torque = custom_tn_df['torque'].max()
                        manual_max_power = custom_tn_df['power_kw'].max()
                        base_speed_idx = custom_tn_df['power_kw'].idxmax()
                        base_speed_calc = custom_tn_df.loc[base_speed_idx, 'rpm']
                        st.info(f"📊 目前載入 [**{t_col}**]：\n最大扭矩 **{manual_peak_torque:.1f} Nm**, 最高轉速 **{manual_max_rpm:.0f} rpm**, 最大功率 **{manual_max_power:.2f} kW**")
                except Exception as e:
                    st.error(f"解析檔案失敗: {e}")

        motor_eff = st.number_input("馬達效率 (%)", min_value=0.0, max_value=100.0, value=90.0, step=1.0)
        
        # ---------- 效率地圖上傳（支援寬格式，靜默處理）----------
        st.markdown("#### 📊 馬達效率地圖 (選用)")
        eff_file = st.file_uploader("上傳效率地圖 CSV", type=["csv"], key="eff_map")
        eff_interpolator = None
        df_eff_converted = None

        if eff_file is not None and SCIPY_AVAILABLE:
            try:
                df_eff_raw = pd.read_csv(eff_file)
                first_col = df_eff_raw.columns[0]
                other_cols = df_eff_raw.columns[1:]
                is_wide_format = False
                if any(key in first_col.lower() for key in ['rpm', '轉速', 'speed']):
                    if any(('%' in col) or ('效率' in col) for col in other_cols):
                        is_wide_format = True
                
                if is_wide_format:
                    if est_mode == '手動輸入' and 'manual_peak_torque' in locals() and manual_peak_torque is not None:
                        default_peak_torque = manual_peak_torque
                    else:
                        default_peak_torque = 18.0
                    peak_torque_input = st.number_input("馬達峰值扭矩 (Nm) - 用於效率地圖換算", 
                                                        min_value=0.1, value=default_peak_torque, step=0.5, key="eff_map_peak_torque")
                    rpm_col = first_col
                    df_eff_raw[rpm_col] = pd.to_numeric(df_eff_raw[rpm_col], errors='coerce')
                    df_eff_raw = df_eff_raw.dropna(subset=[rpm_col])
                    records = []
                    for col in other_cols:
                        match = re.search(r'(\d+(?:\.\d+)?)', col)
                        if match:
                            percent = float(match.group(1))
                        else:
                            continue
                        torque_actual = peak_torque_input * (percent / 100.0)
                        eff_series = df_eff_raw[col].astype(str).str.replace(r'&nbsp;', '', regex=True)
                        eff_series = pd.to_numeric(eff_series, errors='coerce')
                        for idx, row in df_eff_raw.iterrows():
                            rpm_val = row[rpm_col]
                            eff_val = eff_series.iloc[idx]
                            if not np.isnan(rpm_val) and not np.isnan(eff_val) and eff_val >= 0:
                                records.append({'rpm': rpm_val, 'torque': torque_actual, 'efficiency': eff_val})
                    df_eff_converted = pd.DataFrame(records)
                    st.success(f"✅ 成功讀取寬格式效率地圖，共 {len(df_eff_converted)} 筆資料 (扭矩範圍 {df_eff_converted['torque'].min():.1f} ~ {df_eff_converted['torque'].max():.1f} Nm)")
                else:
                    rpm_candidates = [c for c in df_eff_raw.columns if 'rpm' in c.lower() or '轉速' in c]
                    torque_candidates = [c for c in df_eff_raw.columns if 'torque' in c.lower() or '扭矩' in c or '扭力' in c]
                    eff_candidates = [c for c in df_eff_raw.columns if 'eff' in c.lower() or '效率' in c]
                    if rpm_candidates and torque_candidates and eff_candidates:
                        rpm_col = rpm_candidates[0]
                        torque_col = torque_candidates[0]
                        eff_col = eff_candidates[0]
                    else:
                        rpm_col, torque_col, eff_col = df_eff_raw.columns[0], df_eff_raw.columns[1], df_eff_raw.columns[2]
                    df_eff_converted = df_eff_raw[[rpm_col, torque_col, eff_col]].copy()
                    df_eff_converted.columns = ['rpm', 'torque', 'efficiency']
                    df_eff_converted = df_eff_converted.dropna()
                    st.success(f"✅ 成功讀取效率地圖，共 {len(df_eff_converted)} 筆資料")
                
                if df_eff_converted is not None and len(df_eff_converted) > 0:
                    eff_interpolator = build_efficiency_interpolator(df_eff_converted)
                    st.session_state.df_eff_converted = df_eff_converted
                else:
                    st.error("❌ 無有效效率數據")
            except Exception as e:
                st.error(f"❌ 讀取效率地圖失敗: {e}")
        elif eff_file is not None and not SCIPY_AVAILABLE:
            st.error("❌ 需要安裝 SciPy 才能使用效率地圖功能。")

    with st.expander("🔹 齒輪 (Gear)", expanded=True):
        gear_option = st.radio("減速比", ['自動估算', '手動輸入'], index=1)
        if gear_option == '手動輸入':
            gear_ratio = st.number_input("請輸入減速比", min_value=1.0, value=8.7, step=0.5)
        else:
            gear_ratio = None
        gear_eff = st.number_input("齒輪箱效率 (%)", min_value=0.0, max_value=100.0, value=95.0, step=1.0)

    with st.expander("🔹 電池 (Battery)", expanded=True):
        user_battery_energy_kwh = st.number_input("電池總能量 (kWh)", min_value=0.5, value=5.0, step=0.5)
        battery_soc = st.number_input("電池可用 SOC (%)", min_value=0.0, max_value=100.0, value=90.0, step=5.0)

    # 行駛工況設定
    st.header("📁 行駛工況設定")
    st.markdown("上傳 **行駛工況 CSV**（需含時間(s)、車速(km/h)；若無加速度欄位，系統將自動計算）")
    wltc_file = st.file_uploader("選擇行駛工況 CSV 檔案", type=["csv"], key="wltc")
    if wltc_file is not None:
        try:
            df_wltc = pd.read_csv(wltc_file)
            st.success(f"成功讀取行駛工況，共 {len(df_wltc)} 筆資料")
            time_candidates = [c for c in df_wltc.columns if 'time' in c.lower() or 't' in c.lower()]
            speed_candidates = [c for c in df_wltc.columns if 'speed' in c.lower() or 'velocity' in c.lower() or 'v' in c.lower()]
            time_col = st.selectbox("時間欄位 (秒)", df_wltc.columns, index=df_wltc.columns.get_loc(time_candidates[0]) if time_candidates else 0)
            speed_col = st.selectbox("車速欄位 (km/h)", df_wltc.columns, index=df_wltc.columns.get_loc(speed_candidates[0]) if speed_candidates else 1 if len(df_wltc.columns) > 1 else 0)
            df_wltc_clean = df_wltc[[time_col, speed_col]].copy()
            df_wltc_clean.columns = ['time', 'speed_kmh']
            df_wltc_clean = df_wltc_clean.dropna()
            accel_col = next((col for col in df_wltc.columns if 'accel' in col.lower() or 'a' in col.lower()), None)
            if accel_col is not None:
                df_wltc_clean['accel_ms2'] = df_wltc[accel_col]
            else:
                speeds_ms = df_wltc_clean['speed_kmh'].values / 3.6
                times = df_wltc_clean['time'].values
                df_wltc_clean['accel_ms2'] = np.gradient(speeds_ms, times)
            st.session_state.df_wltc_clean = df_wltc_clean
            if 'gear_ratio' in locals() and gear_ratio is not None:
                gear_ratio_val = gear_ratio
            else:
                gear_ratio_val = estimate_gearbox(speed_ms, wheel_radius_m)
            gear_eff_val = gear_eff / 100.0
            df_op = compute_motor_operating_points_from_wltc(
                df_wltc_clean, total_mass, area, cd, fr, wheel_radius_m, gear_ratio_val, gear_eff_val
            )
            st.session_state.df_motor_operating_points = df_op
            st.success("已計算工作點，將在圖1和圖2中疊加顯示。")
        except Exception as e:
            st.error(f"讀取檔案失敗: {e}")
            if "df_wltc_clean" in st.session_state:
                del st.session_state.df_wltc_clean
                del st.session_state.df_motor_operating_points
    else:
        if "df_wltc_clean" in st.session_state:
            del st.session_state.df_wltc_clean
            del st.session_state.df_motor_operating_points

    st.markdown("---")
    st.caption("修改參數後，下方結果會自動更新")


# ================== 計算核心 ==================
if 'gear_ratio' not in locals() or gear_ratio is None:
    gear_ratio = estimate_gearbox(speed_ms, wheel_radius_m)
if 'gear_eff' not in locals():
    gear_eff = 95.0

if voltage is None:
    voltage = 48 if manual_max_power < 20 else 96

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
    if est_mode == '讀取馬達TN曲線' and custom_tn_df is not None:
         base_speed = custom_tn_df.loc[custom_tn_df['power_kw'].idxmax(), 'rpm']
         motor_spec['基速 (rpm)'] = base_speed

rated_power = max_power_kw_used / 2

if desired_range:
    time_h = desired_range / (speed_ms * 0.7 * 3.6)
    battery_spec = estimate_battery(rated_power * 0.7, voltage, duration_h=time_h)
else:
    battery_spec = estimate_battery(rated_power, voltage, duration_h=1.0)

controller_spec = estimate_controller(max_power_kw_used, voltage)
gearbox_spec = {'類型': '固定減速比齒輪箱', '減速比': round(gear_ratio, 2), '效率 (%)': gear_eff}

F_roll_start = total_mass * G * fr
T_motor_start_full = (F_roll_start + total_mass * avg_accel_full) * wheel_radius_m / (gear_ratio * ETA_DRIVE)
T_motor_start_50 = (F_roll_start + total_mass * avg_accel_50) * wheel_radius_m / (gear_ratio * ETA_DRIVE)

n = np.linspace(0, n_max_motor * 1.1, 500)
T_motor_max = np.zeros_like(n)
P_motor_out = np.zeros_like(n)

if est_mode == '讀取馬達TN曲線' and custom_tn_df is not None:
    T_motor_max = np.interp(n, custom_tn_df['rpm'].values, custom_tn_df['torque'].values, right=0)
    P_motor_out = T_motor_max * n / 9550.0
else:
    const_idx = n <= base_speed
    T_motor_max[const_idx] = T_peak
    P_motor_out[const_idx] = T_peak * n[const_idx] / 9550
    power_idx = (n > base_speed) & (n <= n_max_motor)
    T_motor_max[power_idx] = (max_power_kw_used * 1000) / (2 * math.pi * n[power_idx] / 60)
    P_motor_out[power_idx] = max_power_kw_used
    over_idx = n > n_max_motor
    T_motor_max[over_idx] = 0
    P_motor_out[over_idx] = 0

v_from_n = n / gear_ratio * (2 * math.pi * wheel_radius_m) * 3.6 / 60
T_wheel_max = T_motor_max * gear_ratio * ETA_DRIVE
F_wheel_max = T_wheel_max / wheel_radius_m  # 車輪推力 (N)

v_max_motor = n_max_motor / gear_ratio * (2 * math.pi * wheel_radius_m) * 3.6 / 60

motor_rpm_flat, torque_flat, speed_kmh_flat, force_flat = calculate_load_curve(
    total_mass, area, cd, fr, wheel_radius_m, gear_ratio, speed_ms, grade_percent=0, extend_to_vmax=v_max_motor
)
if grade_percent > 0:
    motor_rpm_climb, torque_climb, speed_kmh_climb, force_climb = calculate_load_curve(
        total_mass, area, cd, fr, wheel_radius_m, gear_ratio, speed_ms, grade_percent, extend_to_vmax=v_max_motor
    )
else:
    motor_rpm_climb, torque_climb = None, None

# 車輪推力曲線
F_wheel_flat = force_flat  # 原本 force_flat 已經是總推力 (N)
F_wheel_climb = force_climb if grade_percent > 0 else None

time_acc, speed_acc, disp_acc = simulate_acceleration(
    total_mass, area, cd, fr, wheel_radius_m, gear_ratio, motor_spec, base_speed, T_peak, speed_ms, dt=0.1, custom_tn_df=custom_tn_df
)
actual_0to50 = time_acc[np.argmax(speed_acc >= 50)] if np.any(speed_acc >= 50) else np.inf
actual_full_time = time_acc[np.argmax(speed_acc >= speed_kmh * 0.99)] if np.any(speed_acc >= speed_kmh * 0.99) else np.inf

if "default_motor_spec" not in st.session_state:
    st.session_state.default_motor_spec = motor_spec.copy()
if "default_battery_spec" not in st.session_state:
    st.session_state.default_battery_spec = battery_spec.copy()

# ================== 顯示區 ==================
st.subheader("📋 規格摘要")
with st.expander("📦 馬達規格", expanded=True):
    st.markdown(render_json_with_diff(motor_spec, st.session_state.default_motor_spec), unsafe_allow_html=True)

with st.expander("🏎️ 動態性能表現", expanded=True):
    st.markdown("**起步性能比較**")
    st.metric("0→最高車速起步所需馬達扭矩", f"{T_motor_start_full:.1f} Nm")
    st.metric("0→50 km/h 起步所需馬達扭矩", f"{T_motor_start_50:.1f} Nm")
    if T_motor_start_full <= T_peak and T_motor_start_50 <= T_peak:
        st.success("✅ 馬達峰值扭矩足夠滿足兩種加速需求")
    else:
        st.error(f"❌ 馬達峰值扭矩不足，需增加 {max(0, T_motor_start_full - T_peak, T_motor_start_50 - T_peak):.1f} Nm")
    st.markdown("---")
    col_a, col_b = st.columns(2)
    with col_a:
        st.metric("目標 0→50 km/h", f"{accel_time_0to50:.1f} s")
        st.markdown(f'<p style="color:{"green" if actual_0to50 <= accel_time_0to50 else "red"};">實際: {actual_0to50:.1f} s</p>', unsafe_allow_html=True)
    with col_b:
        st.metric("目標 0→最高車速", f"{accel_time_full:.1f} s")
        st.markdown(f'<p style="color:{"green" if actual_full_time <= accel_time_full else "red"};">實際: {actual_full_time:.1f} s</p>', unsafe_allow_html=True)

with st.expander("🔋 電池 (估算規格)", expanded=False):
    st.markdown(render_battery_with_diff(battery_spec, st.session_state.default_battery_spec), unsafe_allow_html=True)
with st.expander("🎛️ 控制器", expanded=False):
    st.json(controller_spec)
with st.expander("⚙️ 齒輪箱", expanded=False):
    st.json(gearbox_spec)
with st.expander("🔁 轉換係數", expanded=False):
    st.metric("輪上扭矩 / 馬達扭矩", f"{gear_ratio * ETA_DRIVE:.3f}")
    st.metric("車速 (km/h) / 馬達轉速 (rpm)", f"{(2 * math.pi * wheel_radius_m * 60) / (gear_ratio * 1000) * 3.6:.6f}")

st.markdown("---")

# ================== 圖1：馬達 TN (含功率&工作點) ==================
st.markdown("## 📈 圖1：馬達 TN (含功率&工作點)")

x_upper = n_max_motor * 1.1
grid_step = T_peak / 4.0 if T_peak > 0 else 10
y_min_raw = min(0, T_motor_max.min(), torque_flat.min())
if "df_motor_operating_points" in st.session_state:
    min_op = st.session_state.df_motor_operating_points['motor_torque_Nm'].min()
    if not np.isnan(min_op): y_min_raw = min(y_min_raw, min_op)
y_min_torque = math.floor(y_min_raw / grid_step) * grid_step
y_max_torque = T_peak + grid_step
if "df_motor_operating_points" in st.session_state:
    max_op = st.session_state.df_motor_operating_points['motor_torque_Nm'].max()
    if not np.isnan(max_op) and max_op > y_max_torque: y_max_torque = math.ceil(max_op / grid_step) * grid_step

ratio = max_power_kw_used / T_peak if T_peak > 0 else 1
p_min = y_min_torque * ratio
p_max = y_max_torque * ratio

num_ticks = int(round((y_max_torque - y_min_torque) / grid_step)) + 1
y_ticks = [round(v, 2) for v in [y_min_torque + i * grid_step for i in range(num_ticks)] if abs(v - T_peak) > (grid_step*0.1)]
p_ticks = [round(v * ratio, 2) for v in y_ticks]
x_ticks = sorted(list(set([round(v, -1) for v in list(np.linspace(0, x_upper, 6)) + [base_speed, n_max_motor]])))

fig1 = make_subplots(specs=[[{"secondary_y": True}]])
fig1.add_trace(go.Scatter(x=n, y=T_motor_max, mode='lines', name='馬達最大扭矩', line=dict(color='dodgerblue', width=3)), secondary_y=False)
fig1.add_trace(go.Scatter(x=motor_rpm_flat, y=torque_flat, mode='lines', name='平路負載線', line=dict(color='red', width=3, dash='dash')), secondary_y=False)
if motor_rpm_climb is not None:
    fig1.add_trace(go.Scatter(x=motor_rpm_climb, y=torque_climb, mode='lines', name=f'爬坡負載線', line=dict(color='green', width=3, dash='dot')), secondary_y=False)
fig1.add_trace(go.Scatter(x=n, y=P_motor_out, mode='lines', name='馬達功率', line=dict(color='gold', width=2, dash='solid')), secondary_y=True)

fig1.add_trace(go.Scatter(x=[0, base_speed], y=[T_peak, T_peak], mode='markers', name='關鍵點', marker=dict(color='dodgerblue', size=10)), secondary_y=False)
T_at_max_n = (max_power_kw_used * 1000) / (2 * math.pi * n_max_motor / 60) if n_max_motor > 0 else 0
fig1.add_trace(go.Scatter(x=[n_max_motor], y=[T_at_max_n], mode='markers', name='極速點', marker=dict(color='purple', size=10)), secondary_y=False)

design_rpm = speed_ms * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio
fig1.add_vline(x=design_rpm, line_width=2, line_dash="dash", line_color="orange", opacity=0.9)
T_at_design = np.interp(design_rpm, n, T_motor_max) if design_rpm <= n_max_motor else 0
fig1.add_trace(go.Scatter(x=[design_rpm], y=[T_at_design], mode='markers', name='目標轉速', marker=dict(color='orange', size=10)), secondary_y=False)

if "df_motor_operating_points" in st.session_state:
    df_op = st.session_state.df_motor_operating_points
    fig1.add_trace(go.Scatter(x=df_op['motor_rpm'], y=df_op['motor_torque_Nm'], mode='markers', marker=dict(size=4, color='cyan', opacity=0.6, symbol='circle'), name='工作點', showlegend=True), secondary_y=False)

# ---------- 平路負載線與馬達最大扭矩曲線的交點（含數值標示）----------
intersections_flat = find_intersection(n, T_motor_max, motor_rpm_flat, torque_flat)
for idx, (x_cross, y_cross) in enumerate(intersections_flat):
    # 標記點
    fig1.add_trace(go.Scatter(
        x=[x_cross], y=[y_cross], mode='markers',
        marker=dict(color='red', size=10, symbol='x'),
        name='平路交點', showlegend=False,
        text=f"平路交點<br>轉速: {x_cross:.0f} rpm<br>扭矩: {y_cross:.1f} Nm",
        hoverinfo='text'
    ))
    # 數值標示（偏移避免重疊）
    fig1.add_annotation(
        x=x_cross, y=y_cross,
        text=f"<b>{x_cross:.0f} rpm, {y_cross:.1f} Nm</b>",
        showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=1,
        ax=40, ay=-30 - idx*20,
        font=dict(size=10, color='red'),
        bgcolor="rgba(0,0,0,0.6)", bordercolor='red', borderwidth=1
    )

# ---------- 爬坡負載線與馬達最大扭矩曲線的交點（含數值標示）----------
if motor_rpm_climb is not None:
    intersections_climb = find_intersection(n, T_motor_max, motor_rpm_climb, torque_climb)
    for idx, (x_cross, y_cross) in enumerate(intersections_climb):
        fig1.add_trace(go.Scatter(
            x=[x_cross], y=[y_cross], mode='markers',
            marker=dict(color='green', size=10, symbol='x'),
            name='爬坡交點', showlegend=False,
            text=f"爬坡交點<br>轉速: {x_cross:.0f} rpm<br>扭矩: {y_cross:.1f} Nm",
            hoverinfo='text'
        ))
        fig1.add_annotation(
            x=x_cross, y=y_cross,
            text=f"<b>{x_cross:.0f} rpm, {y_cross:.1f} Nm</b>",
            showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=1,
            ax=-60, ay=-30 - idx*20,
            font=dict(size=10, color='green'),
            bgcolor="rgba(0,0,0,0.6)", bordercolor='green', borderwidth=1
        )

fig1.update_yaxes(title_text="扭矩 (Nm)", secondary_y=False, range=[y_min_torque, y_max_torque], tickvals=y_ticks, tickfont=dict(color='white'), zeroline=True, zerolinecolor='gray', zerolinewidth=1.5)
fig1.update_yaxes(title_text="功率 (kW)", secondary_y=True, range=[p_min, p_max], tickvals=p_ticks, tickfont=dict(color='white'), zeroline=True, zerolinecolor='gray', zerolinewidth=1.5, showgrid=False)
fig1.update_xaxes(title_text="轉速 (rpm)", range=[0, x_upper], tickvals=x_ticks, tickfont=dict(color='white'), zeroline=True, zerolinecolor='gray', zerolinewidth=1.5)

fig1.add_annotation(x=0, y=T_peak, xref="paper", yref="y", text=f"<b>{T_peak:.1f}</b>", showarrow=False, xanchor="right", xshift=-15, font=dict(color="dodgerblue", size=14), bgcolor="rgba(26,28,35,0.9)", bordercolor="dodgerblue", borderwidth=1, borderpad=4)
fig1.add_annotation(x=1, y=max_power_kw_used, xref="paper", yref="y2", text=f"<b>{max_power_kw_used:.2f}</b>", showarrow=False, xanchor="left", xshift=15, font=dict(color="gold", size=14), bgcolor="rgba(26,28,35,0.9)", bordercolor="gold", borderwidth=1, borderpad=4)
fig1.add_annotation(x=base_speed, y=T_peak, xref="x", yref="y", text=f"<b>基速: {base_speed:.0f} rpm</b>", showarrow=True, arrowhead=2, arrowcolor="green", arrowsize=1, arrowwidth=2, ax=0, ay=-45, font=dict(color="lightgreen", size=12), bgcolor="rgba(26,28,35,0.9)", bordercolor="green", borderwidth=1, borderpad=3)
fig1.add_annotation(x=design_rpm, y=T_at_design, xref="x", yref="y", text=f"<b>目標: {design_rpm:.0f} rpm</b>", showarrow=True, arrowhead=2, arrowcolor="orange", arrowsize=1, arrowwidth=2, ax=-50, ay=-70, font=dict(color="orange", size=12), bgcolor="rgba(26,28,35,0.9)", bordercolor="orange", borderwidth=1, borderpad=3)
fig1.add_annotation(x=n_max_motor, y=T_at_max_n, xref="x", yref="y", text=f"<b>極速: {n_max_motor:.0f} rpm<br>{T_at_max_n:.1f} Nm</b>", showarrow=True, arrowhead=2, arrowcolor="purple", arrowsize=1, arrowwidth=2, ax=60, ay=-45, font=dict(color="#d8b4e2", size=12), bgcolor="rgba(26,28,35,0.9)", bordercolor="purple", borderwidth=1, borderpad=3)

fig1.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5), margin=dict(l=80, r=110, t=100, b=20), height=550)
st.plotly_chart(fig1, use_container_width=True)
st.markdown("---")


# ================== 圖2：車輪推力(含工作點) ==================
st.markdown("## 📈 圖2：車輪推力(含工作點)")

idx_design = np.argmin(np.abs(speed_kmh_flat - speed_kmh))
F_design_flat = F_wheel_flat[idx_design]
F_at_vmax = np.interp(v_max_motor, v_from_n, F_wheel_max) if v_max_motor <= v_from_n.max() else 0
F_wheel_peak = F_wheel_max.max()

grid_step_force = F_wheel_peak / 4.0 if F_wheel_peak > 0 else 10
y_min_raw_f = min(0, F_wheel_max.min(), F_wheel_flat.min())
if "df_motor_operating_points" in st.session_state:
    wheel_force_op = st.session_state.df_motor_operating_points['wheel_torque_Nm'] / wheel_radius_m
    min_op_f = wheel_force_op.min()
    if not np.isnan(min_op_f): y_min_raw_f = min(y_min_raw_f, min_op_f)
y_min_force = math.floor(y_min_raw_f / grid_step_force) * grid_step_force

y_max_force = F_wheel_peak + grid_step_force
if "df_motor_operating_points" in st.session_state:
    max_op_f = wheel_force_op.max()
    if not np.isnan(max_op_f) and max_op_f > y_max_force: y_max_force = math.ceil(max_op_f / grid_step_force) * grid_step_force

num_ticks_f = int(round((y_max_force - y_min_force) / grid_step_force)) + 1
y_ticks_f = [round(v, 2) for v in [y_min_force + i * grid_step_force for i in range(num_ticks_f)] if abs(v - F_wheel_peak) > (grid_step_force*0.1)]

fig2 = go.Figure()
fig2.add_trace(go.Scatter(x=v_from_n, y=F_wheel_max, mode='lines', name='最大車輪推力', line=dict(color='dodgerblue', width=3)))
fig2.add_trace(go.Scatter(x=speed_kmh_flat, y=F_wheel_flat, mode='lines', name='平路負載線', line=dict(color='red', width=3, dash='dash')))
if F_wheel_climb is not None:
    fig2.add_trace(go.Scatter(x=speed_kmh_climb, y=F_wheel_climb, mode='lines', name=f'爬坡負載線', line=dict(color='green', width=3, dash='dot')))

fig2.add_vline(x=speed_kmh, line_width=2, line_dash="dash", line_color="orange", opacity=0.9)
fig2.add_trace(go.Scatter(x=[speed_kmh, v_max_motor], y=[F_design_flat, F_at_vmax], mode='markers', name='關鍵點', marker=dict(color='orange', size=10)))

if "df_motor_operating_points" in st.session_state:
    df_op = st.session_state.df_motor_operating_points
    wheel_force_op = df_op['wheel_torque_Nm'] / wheel_radius_m
    fig2.add_trace(go.Scatter(x=df_op['speed_kmh'], y=wheel_force_op, mode='markers', marker=dict(size=4, color='cyan', opacity=0.6, symbol='circle'), name='工作點', showlegend=True))

# ---------- 平路負載線與最大車輪推力曲線的交點（含數值標示）----------
intersections_flat_force = find_intersection(v_from_n, F_wheel_max, speed_kmh_flat, F_wheel_flat)
for idx, (x_cross, y_cross) in enumerate(intersections_flat_force):
    fig2.add_trace(go.Scatter(
        x=[x_cross], y=[y_cross], mode='markers',
        marker=dict(color='red', size=10, symbol='x'),
        name='平路交點', showlegend=False,
        text=f"平路交點<br>車速: {x_cross:.1f} km/h<br>推力: {y_cross:.0f} N",
        hoverinfo='text'
    ))
    fig2.add_annotation(
        x=x_cross, y=y_cross,
        text=f"<b>{x_cross:.1f} km/h, {y_cross:.0f} N</b>",
        showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=1,
        ax=40, ay=-30 - idx*20,
        font=dict(size=10, color='red'),
        bgcolor="rgba(0,0,0,0.6)", bordercolor='red', borderwidth=1
    )

# ---------- 爬坡負載線與最大車輪推力曲線的交點（含數值標示）----------
if F_wheel_climb is not None:
    intersections_climb_force = find_intersection(v_from_n, F_wheel_max, speed_kmh_climb, F_wheel_climb)
    for idx, (x_cross, y_cross) in enumerate(intersections_climb_force):
        fig2.add_trace(go.Scatter(
            x=[x_cross], y=[y_cross], mode='markers',
            marker=dict(color='green', size=10, symbol='x'),
            name='爬坡交點', showlegend=False,
            text=f"爬坡交點<br>車速: {x_cross:.1f} km/h<br>推力: {y_cross:.0f} N",
            hoverinfo='text'
        ))
        fig2.add_annotation(
            x=x_cross, y=y_cross,
            text=f"<b>{x_cross:.1f} km/h, {y_cross:.0f} N</b>",
            showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=1,
            ax=-60, ay=-30 - idx*20,
            font=dict(size=10, color='green'),
            bgcolor="rgba(0,0,0,0.6)", bordercolor='green', borderwidth=1
        )

x_max = max(v_max_motor, speed_kmh) * 1.15 if max(v_max_motor, speed_kmh) > 0 else 100
fig2.update_yaxes(title_text="車輪推力 (N)", range=[y_min_force, y_max_force], tickvals=y_ticks_f, tickfont=dict(color='white'), zeroline=True, zerolinecolor='gray', zerolinewidth=1.5)
fig2.update_xaxes(title_text="車速 (km/h)", range=[0, x_max], tickfont=dict(color='white'), zeroline=True, zerolinecolor='gray', zerolinewidth=1.5)

fig2.add_annotation(x=0, y=F_wheel_peak, xref="x", yref="y", text=f"<b>{F_wheel_peak:.0f}</b>", showarrow=False, xanchor="right", xshift=-15, font=dict(color="dodgerblue", size=14), bgcolor="rgba(26,28,35,0.9)", bordercolor="dodgerblue", borderwidth=1, borderpad=4)
fig2.add_annotation(x=speed_kmh, y=F_design_flat, xref="x", yref="y", text=f"<b>目標: {speed_kmh:.0f} km/h</b>", showarrow=True, arrowhead=2, arrowcolor="orange", arrowsize=1, arrowwidth=2, ax=-55, ay=-65, font=dict(color="orange", size=12), bgcolor="rgba(26,28,35,0.9)", bordercolor="orange", borderwidth=1, borderpad=3)
fig2.add_annotation(x=v_max_motor, y=F_at_vmax, xref="x", yref="y", text=f"<b>極速: {v_max_motor:.0f} km/h<br>{F_at_vmax:.0f} N</b>", showarrow=True, arrowhead=2, arrowcolor="purple", arrowsize=1, arrowwidth=2, ax=65, ay=-45, font=dict(color="#d8b4e2", size=12), bgcolor="rgba(26,28,35,0.9)", bordercolor="purple", borderwidth=1, borderpad=3)

fig2.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5), margin=dict(l=80, r=110, t=100, b=20), height=550)
st.plotly_chart(fig2, use_container_width=True)
st.markdown("---")


# ================== 圖3：加速性能（速度與位移 vs 時間）==================
st.markdown("## 📈 圖3：加速性能（速度與位移 vs 時間）")
st.caption("藍色實線為車速隨時間變化，紅色虛線為位移隨時間變化。垂直線標註實際達到50 km/h和最高車速的時間，以及目標加速時間。")

fig3 = make_subplots(specs=[[{"secondary_y": True}]])
fig3.add_trace(go.Scatter(x=time_acc, y=speed_acc, mode='lines', name='車速 (km/h)', line=dict(color='dodgerblue', width=3)), secondary_y=False)
fig3.add_trace(go.Scatter(x=time_acc, y=disp_acc, mode='lines', name='位移 (m)', line=dict(color='red', width=2, dash='dash')), secondary_y=True)

idx_50 = np.argmax(speed_acc >= 50)
if idx_50 > 0:
    t_50 = time_acc[idx_50]
    fig3.add_vline(x=t_50, line_width=1, line_dash="dot", line_color="orange", opacity=0.7)
    fig3.add_annotation(x=t_50, y=50, text=f"實際50km/h @ {t_50:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=-30)
idx_max = np.argmax(speed_acc >= speed_kmh * 0.99)
if idx_max > 0:
    t_max = time_acc[idx_max]
    fig3.add_vline(x=t_max, line_width=1, line_dash="dot", line_color="green", opacity=0.7)
    fig3.add_annotation(x=t_max, y=speed_kmh, text=f"實際{int(speed_kmh)}km/h @ {t_max:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=30)
fig3.add_vline(x=accel_time_0to50, line_width=1, line_dash="dot", line_color="purple", opacity=0.7)
fig3.add_annotation(x=accel_time_0to50, y=50, text=f"目標50km/h @ {accel_time_0to50:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=-50, font=dict(color="purple"))
fig3.add_vline(x=accel_time_full, line_width=1, line_dash="dot", line_color="brown", opacity=0.7)
fig3.add_annotation(x=accel_time_full, y=speed_kmh, text=f"目標最高車速 @ {accel_time_full:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=-70, font=dict(color="brown"))

fig3.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5), margin=dict(l=20, r=20, t=40, b=20), height=400)
fig3.update_xaxes(title_text="時間 (秒)")
fig3.update_yaxes(title_text="車速 (km/h)", secondary_y=False)
fig3.update_yaxes(title_text="位移 (m)", secondary_y=True)
st.plotly_chart(fig3, use_container_width=True)
st.caption("💡 提示：圖中紫色虛線為目標 0→50 km/h 加速時間，棕色虛線為目標 0→最高車速加速時間。")
st.markdown("---")

# ================== 圖4 區塊 ==================
if "df_wltc_clean" in st.session_state:
    st.markdown("## 📈 圖4：行駛工況（車速與加速度 vs 時間）")
    df_wltc_plot = st.session_state.df_wltc_clean
    fig4 = make_subplots(specs=[[{"secondary_y": True}]])
    fig4.add_trace(go.Scatter(x=df_wltc_plot['time'], y=df_wltc_plot['speed_kmh'], mode='lines', name='車速', line=dict(color='dodgerblue', width=2)), secondary_y=False)
    fig4.add_trace(go.Scatter(x=df_wltc_plot['time'], y=df_wltc_plot['accel_ms2'], mode='lines', name='加速度', line=dict(color='red', width=2, dash='dash')), secondary_y=True)
    fig4.update_layout(height=400, margin=dict(l=20, r=20, t=40, b=20))
    st.plotly_chart(fig4, use_container_width=True)
    st.markdown("---")

# ================== 圖5：電池能耗分析與里程預估 ==================
if "df_wltc_clean" in st.session_state and SCIPY_AVAILABLE:
    st.markdown("## 📈 圖5：電池能耗分析與里程預估")
    st.caption("本模組基於行駛工況積分計算能耗，並根據側邊欄設定的**電池容量與 SOC** 進行續航里程預估。")
    
    df_energy = st.session_state.df_wltc_clean.copy()
    df_operating = st.session_state.df_motor_operating_points if "df_motor_operating_points" in st.session_state else None
    
    gear_eff_val = gear_eff if 'gear_eff' in locals() else 95.0
    motor_eff_val = motor_eff if 'motor_eff' in locals() else 90.0
    
    wh_per_km_batt, wh_per_km_wheel, total_dist_km, total_energy_wh, drive_energy_wh, regen_energy_wh, times, power_batt_w, energy_cum = compute_theoretical_energy_consumption(
        df_energy, total_mass, area, cd, fr, gear_eff_val, motor_eff_val, eff_interpolator, df_operating
    )
    
    usable_energy_wh = user_battery_energy_kwh * 1000 * (battery_soc / 100.0)
    estimated_range_km = usable_energy_wh / wh_per_km_batt if wh_per_km_batt > 0 else 0
    
    c1, c2, c3 = st.columns(3)
    c1.metric("🏁 工況總里程", f"{total_dist_km:.3f} km")
    c2.metric("⚙️ 輪上能耗", f"{wh_per_km_wheel:.2f} Wh/km")
    c3.metric("🔋 電池端輸出能耗", f"{wh_per_km_batt:.2f} Wh/km")
    
    st.markdown("#### 🔋 續航里程估算結果")
    cc1, cc2 = st.columns(2)
    cc1.metric("⚡ 可用總能量", f"{usable_energy_wh:.1f} Wh")
    cc2.metric("🎯 理論預估里程", f"{estimated_range_km:.1f} km", delta=f"基於 {user_battery_energy_kwh:.1f}kWh 電池")
    
    fig5 = make_subplots(specs=[[{"secondary_y": True}]])
    power_pos = np.maximum(power_batt_w, 0)
    power_neg = np.minimum(power_batt_w, 0)
    fig5.add_trace(go.Scatter(x=times, y=power_pos, mode='lines', fill='tozeroy', name='驅動輸出 (W)', line=dict(color='crimson', width=1)), secondary_y=False)
    fig5.add_trace(go.Scatter(x=times, y=power_neg, mode='lines', fill='tozeroy', name='動能回收 (W)', line=dict(color='seagreen', width=1)), secondary_y=False)
    fig5.add_trace(go.Scatter(x=times, y=energy_cum, mode='lines', name='累積淨耗電 (Wh)', line=dict(color='gold', width=3, dash='solid')), secondary_y=True)
    fig5.update_xaxes(title_text="時間 (秒)", zeroline=True, zerolinecolor='gray')
    fig5.update_yaxes(title_text="瞬時電池功率 (W)", secondary_y=False, zeroline=True, zerolinecolor='gray')
    fig5.update_yaxes(title_text="累積耗電 (Wh)", secondary_y=True, zeroline=False)
    fig5.update_layout(height=450, margin=dict(l=20, r=20, t=40, b=20), legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5))
    st.plotly_chart(fig5, use_container_width=True)
    
    with st.expander("📐 理論能耗計算公式說明", expanded=True):
        st.markdown(r"""
        **1. 行駛阻力 (N)**  
        $$
        F_{\text{total}} = F_{\text{roll}} + F_{\text{aero}} + F_{\text{acc}}
        $$  
        - 滾動阻力：$F_{\text{roll}} = m g f_r$  
        - 空氣阻力：$F_{\text{aero}} = \frac{1}{2} \rho C_d A v^2$  
        - 加速阻力：$F_{\text{acc}} = m a$  

        **2. 輪上功率 (W)**  
        $$
        P_{\text{wheel}} = F_{\text{total}} \cdot v
        $$

        **3. 馬達效率查表 (效率地圖)**  
        根據馬達轉速 $n$ 與扭力 $T$，利用二維線性插值取得即時效率 $\eta_{\text{motor}}(n, T)$。  
        若未上傳效率地圖，則使用固定效率 $\eta_{\text{motor}}$。

        **4. 總效率 (驅動/回收)**  
        - 驅動模式：$\eta_{\text{total}} = \eta_{\text{gear}} \times \eta_{\text{motor}}$  
        - 回收模式：$\eta_{\text{total}} = \eta_{\text{gear}} \times \eta_{\text{motor}}$（回收時功率反向，效率仍相乘）

        **5. 電池端功率 (W)**  
        - 驅動時：$P_{\text{batt}} = P_{\text{wheel}} / \eta_{\text{total}}$  
        - 回收時：$P_{\text{batt}} = P_{\text{wheel}} \times \eta_{\text{total}}$（$P_{\text{wheel}}<0$ 表示減速回收）

        **6. 累積能耗與里程估算**  
        - 累積耗電 (Wh)：$E_{\text{batt}} = \frac{1}{3600} \int P_{\text{batt}} \, dt$  
        - 行駛距離 (km)：$D = \frac{1}{1000} \int v \, dt$  
        - 每公里能耗 (Wh/km)：$\text{EC} = \frac{E_{\text{batt}}}{D}$  
        - 預估里程 (km)：$\text{Range} = \frac{E_{\text{usable}}}{\text{EC}}$，其中 $E_{\text{usable}} = E_{\text{battery}} \times \text{SOC} / 100$

        **7. 離散積分方法**  
        使用 SciPy 的 `integrate.trapezoid` 或 NumPy 的 `np.trapz` 進行梯形法數值積分。
        """)

        st.markdown("---")
        st.markdown("### 📊 本次工況實際計算範例")
        
        try:
            m_kg = total_mass
            fr_val = fr
            cd_val = cd
            area_val = area
            gear_ratio_val = gear_ratio
            gear_eff_val = gear_eff / 100.0
            battery_kwh = user_battery_energy_kwh
            soc = battery_soc
            
            dist_km = total_dist_km
            wh_km_batt = wh_per_km_batt
            drive_wh = drive_energy_wh
            regen_wh = regen_energy_wh
            net_wh = total_energy_wh
            range_km = estimated_range_km
            
            peak_torque_val = None
            if 'peak_torque_input' in locals() and peak_torque_input is not None:
                peak_torque_val = peak_torque_input
            elif 'manual_peak_torque' in locals() and manual_peak_torque is not None:
                peak_torque_val = manual_peak_torque
            else:
                peak_torque_val = T_peak if 'T_peak' in locals() else 0
            
            use_eff_map = (eff_interpolator is not None and df_eff_converted is not None)
            
            st.markdown(f"""
            **🔧 參數盤點與假設**  
            - **總質量**: {m_kg} kg (車重 {weight} kg + 載重 {load} kg)  
            - **阻力設定**: $C_d = {cd_val}$, 迎風面積 $A = {area_val} \, \\text{{m}}^2$, 滾阻 $f_r = {fr_val}$  
            - **傳動幾何**: 輪胎半徑 {wheel_radius_m:.4f} m, 齒輪比 {gear_ratio_val:.1f}, 齒輪效率 {gear_eff}%  
            """)
            
            if use_eff_map:
                if 'df_eff_converted' in st.session_state:
                    torque_min = st.session_state.df_eff_converted['torque'].min()
                    torque_max = st.session_state.df_eff_converted['torque'].max()
                    st.markdown(f"""
            - **馬達效率地圖解析**: 將上傳的效率地圖 CSV 中「轉速」欄位與各扭矩百分比欄位（10%~100%）自動對應到您提供的**峰值扭矩 {peak_torque_val:.1f} Nm** 的百分比（即 Y 軸: {torque_min:.1f} Nm ~ {torque_max:.1f} Nm），建立 2D 效率插值網格。  
                    """)
                else:
                    st.markdown(f"""
            - **馬達效率地圖解析**: 將上傳的效率地圖 CSV 中「轉速」欄位與各扭矩百分比欄位（10%~100%）自動對應到您提供的**峰值扭矩 {peak_torque_val:.1f} Nm** 的百分比，建立 2D 效率插值網格。  
                    """)
            else:
                st.markdown(f"""
            - **馬達效率**: 使用固定效率 {motor_eff_val}% (未上傳效率地圖)。  
                """)
            
            st.markdown(f"""
            - **可用電池能量**: {battery_kwh} kWh × {soc}% = {battery_kwh * soc / 100:.2f} kWh (即 {battery_kwh * soc / 100 * 1000:.0f} Wh)  

            **📈 [行駛工況逐秒積分運算結果]**  
            - **工況單圈總距離**: {dist_km:.3f} km  
            - **電池端輸出能耗**: {wh_km_batt:.2f} Wh/km (對應 {net_wh:.1f} Wh)  
            - **回收充電 (電池端)**: {regen_wh / dist_km:.2f} Wh/km (總計 {regen_wh:.2f} Wh)  
            - **理論預估續航里程**: {range_km:.1f} 公里  

            **✍️ 計算過程舉例（選取某個時間點，例如最高功率點）**  
            """)
            
            if 'times' in locals() and len(times) > 0 and 'power_batt_w' in locals():
                pos_power_mask = power_batt_w > 0
                if np.any(pos_power_mask):
                    idx_max = np.argmax(power_batt_w[pos_power_mask])
                    orig_idx = np.where(pos_power_mask)[0][idx_max]
                    t_ex = times[orig_idx]
                    v_ex = df_energy['speed_kmh'].iloc[orig_idx]
                    a_ex = df_energy['accel_ms2'].iloc[orig_idx]
                    v_ms = v_ex / 3.6
                    F_roll = m_kg * G * fr_val
                    F_aero = 0.5 * RHO * cd_val * area_val * v_ms**2
                    F_acc = m_kg * a_ex
                    F_total = F_roll + F_aero + F_acc
                    P_wheel_ex = F_total * v_ms
                    if df_operating is not None:
                        rpm_ex = df_operating['motor_rpm'].iloc[orig_idx]
                        torque_ex = df_operating['motor_torque_Nm'].iloc[orig_idx]
                    else:
                        rpm_ex = v_ms * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio_val
                        torque_ex = F_total * wheel_radius_m / (gear_ratio_val * gear_eff_val)
                    
                    if use_eff_map and eff_interpolator is not None:
                        eff_motor = eff_interpolator(rpm_ex, torque_ex)[0]
                        st.markdown(f"""
                    在時間 **t = {t_ex:.1f} s** 時，車速 {v_ex:.1f} km/h，加速度 {a_ex:.2f} m/s²。  
                    - 行駛阻力：$F_\\text{{roll}} = {m_kg} \\times 9.81 \\times {fr_val} = {F_roll:.1f}\\,\\text{{N}}$，  
                      $F_\\text{{aero}} = 0.5 \\times 1.2 \\times {cd_val} \\times {area_val} \\times ({v_ms:.2f})^2 = {F_aero:.1f}\\,\\text{{N}}$，  
                      $F_\\text{{acc}} = {m_kg} \\times {a_ex:.3f} = {F_acc:.1f}\\,\\text{{N}}$  
                    - 總阻力：$F_\\text{{total}} = {F_total:.1f}\\,\\text{{N}}$  
                    - 輪上功率：$P_\\text{{wheel}} = {F_total:.1f} \\times {v_ms:.2f} = {P_wheel_ex:.0f}\\,\\text{{W}}$  
                    - 馬達工作點：轉速 {rpm_ex:.0f} rpm，扭矩 {torque_ex:.1f} Nm  
                    - 馬達效率（查表）：$\\eta_\\text{{motor}} = {eff_motor:.1f}\\%$  
                    - 總效率：$\\eta_\\text{{total}} = {gear_eff}% \\times {eff_motor:.1f}\\% = {gear_eff_val * (eff_motor/100) * 100:.1f}\\%$  
                    - 電池端功率：$P_\\text{{batt}} = {P_wheel_ex:.0f} / {gear_eff_val * (eff_motor/100):.3f} = {power_batt_w[orig_idx]:.0f}\\,\\text{{W}}$  
                    - 該點瞬間能耗率即為 $P_\\text{{batt}}$，對全工況時間積分後得到總耗電 {net_wh:.1f} Wh。
                        """)
                    else:
                        st.markdown(f"""
                    在時間 **t = {t_ex:.1f} s** 時，車速 {v_ex:.1f} km/h，加速度 {a_ex:.2f} m/s²。  
                    - 行駛阻力：$F_\\text{{roll}} = {F_roll:.1f}\\,\\text{{N}}$，$F_\\text{{aero}} = {F_aero:.1f}\\,\\text{{N}}$，$F_\\text{{acc}} = {F_acc:.1f}\\,\\text{{N}}$  
                    - 總阻力：$F_\\text{{total}} = {F_total:.1f}\\,\\text{{N}}$  
                    - 輪上功率：$P_\\text{{wheel}} = {P_wheel_ex:.0f}\\,\\text{{W}}$  
                    - 固定馬達效率 {motor_eff_val}%，總效率 $\\eta_\\text{{total}} = {gear_eff}% \\times {motor_eff_val}% = {gear_eff_val * (motor_eff_val/100) * 100:.1f}\\%$  
                    - 電池端功率：$P_\\text{{batt}} = {P_wheel_ex:.0f} / {(gear_eff_val * (motor_eff_val/100)):.3f} = {power_batt_w[orig_idx]:.0f}\\,\\text{{W}}$  
                    - 該點瞬間能耗率即為 $P_\\text{{batt}}$，對全工況時間積分後得到總耗電 {net_wh:.1f} Wh。
                        """)
                else:
                    st.markdown("> 無法找到典型的驅動功率點，跳過詳細計算舉例。")
            else:
                st.markdown("> 由於缺乏逐秒數據，無法展示詳細的時間點計算過程。")
                
        except Exception as e:
            st.markdown(f"> ⚠️ 動態計算範例時發生錯誤：{e}，請確認所有參數已正確設定。")
    
    st.markdown("---")

# ================== 圖6：馬達效率地圖與工作點（獨立顯示）==================
if "df_wltc_clean" in st.session_state and SCIPY_AVAILABLE and eff_interpolator is not None and 'df_eff_converted' in st.session_state:
    df_eff_final = st.session_state.df_eff_converted
    df_operating = st.session_state.df_motor_operating_points if "df_motor_operating_points" in st.session_state else None
    if df_operating is not None:
        st.markdown("## 📈 圖6：馬達效率地圖與工作點分佈")
        st.caption("等高線為效率地圖（%），散點為行駛工況下的馬達操作點，顏色代表該點的效率值。紅色點表示該工作點所需扭矩超過馬達極限（無法滿足）。")
        
        rpm_vals = df_eff_final['rpm'].values
        torque_vals = df_eff_final['torque'].values
        eff_vals = df_eff_final['efficiency'].values
        
        rpm_grid = np.linspace(rpm_vals.min(), rpm_vals.max(), 100)
        torque_grid = np.linspace(torque_vals.min(), torque_vals.max(), 100)
        Rpm_grid, Torque_grid = np.meshgrid(rpm_grid, torque_grid)
        points = np.vstack((rpm_vals, torque_vals)).T
        Eff_grid = griddata(points, eff_vals, (Rpm_grid, Torque_grid), method='cubic')
        
        work_rpm = df_operating['motor_rpm'].values
        work_torque = df_operating['motor_torque_Nm'].values
        eff_at_points = eff_interpolator(work_rpm, work_torque)
        
        max_torque_interp = np.interp(work_rpm, n, T_motor_max, left=0, right=0)
        exceed_mask = work_torque > (max_torque_interp + 1e-3)
        
        fig6 = go.Figure()
        fig6.add_trace(go.Contour(
            x=rpm_grid, y=torque_grid, z=Eff_grid,
            colorscale='Viridis', opacity=0.8,
            contours=dict(coloring='heatmap'),
            colorbar=dict(title="效率 (%)", x=1.02, len=0.8),
            name="效率地圖"
        ))
        
        normal_mask = ~exceed_mask
        if np.any(normal_mask):
            fig6.add_trace(go.Scatter(
                x=work_rpm[normal_mask], y=work_torque[normal_mask],
                mode='markers', marker=dict(size=5, color=eff_at_points[normal_mask], colorscale='Viridis',
                                            colorbar=dict(title="工作點效率 (%)", x=1.15, len=0.8), showscale=True),
                name='工作點 (可滿足)', text=[f"rpm: {r:.0f}<br>Torque: {t:.1f}<br>Eff: {e:.1f}%" 
                                             for r, t, e in zip(work_rpm[normal_mask], work_torque[normal_mask], eff_at_points[normal_mask])],
                hoverinfo='text'
            ))
        
        if np.any(exceed_mask):
            fig6.add_trace(go.Scatter(
                x=work_rpm[exceed_mask], y=work_torque[exceed_mask],
                mode='markers', marker=dict(size=6, color='red', symbol='x'),
                name='工作點 (超出馬達極限)', text=[f"rpm: {r:.0f}<br>Torque: {t:.1f}<br>⚠️ 超出馬達極限 (最大 {max_t:.1f} Nm)" 
                                                   for r, t, max_t in zip(work_rpm[exceed_mask], work_torque[exceed_mask], max_torque_interp[exceed_mask])],
                hoverinfo='text'
            ))
        
        fig6.update_xaxes(title_text="馬達轉速 (rpm)")
        fig6.update_yaxes(title_text="馬達扭矩 (Nm)")
        fig6.update_layout(height=500, margin=dict(l=20, r=140, t=40, b=20),
                          legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5))
        st.plotly_chart(fig6, use_container_width=True)
        st.markdown("---")