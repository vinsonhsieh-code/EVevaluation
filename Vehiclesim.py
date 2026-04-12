import streamlit as st
import pandas as pd
import numpy as np
import math
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from io import BytesIO

# ================== 檢查 SciPy 可用性 ==================
try:
    from scipy import integrate
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    st.warning("⚠️ 未安裝 SciPy，理論能耗計算功能將無法使用。請執行 `pip install scipy` 以啟用此功能。")

# ================== 版本歷史記錄 ==================
st.markdown("""
<details>
<summary>📜 版本歷史記錄</summary>

**v2.2 (2026-04-12) - 穩定除錯與圖形強化版**
- 於圖1及圖2新增負載線與最大扭矩曲線之交點標示，並精確顯示數值。
- 修正 `np.gradient` 導致舊版 SciPy 當機的問題，改用相容性最高的 `np.trapz` 與自製防呆加速度函數。
- 引入強制型別轉換 `pd.to_numeric`，防止上傳的 CSV 夾帶字串導致數學運算崩潰。
- 新增圖5「理論能耗計算」：納入齒輪與馬達效率損耗，區分「輪上能耗」與「電池端理論能耗」。

**v2.1 (2026-04-11)**
- 新增「讀取馬達TN曲線」模式。

**v2.0 (2026-04-06)**
- 整合 WLTC 分析、反向動力學計算。
</details>
""", unsafe_allow_html=True)


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

# ================== 輔助函數 (防彈計算) ==================
def get_safe_accel(speeds_ms, times):
    """ 安全計算加速度，避免 dt=0 或 NaN 導致除以零當機 """
    dt = np.diff(times)
    dt[dt <= 0] = 0.01  # 防禦機制：若時間差小於等於0，強制給一個微小值
    dv = np.diff(speeds_ms)
    accels = np.zeros_like(times)
    accels[1:] = dv / dt
    accels[0] = accels[1] if len(accels) > 1 else 0
    return accels

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

def simulate_acceleration(mass, area, cd, fr, wheel_radius_m, gear_ratio, motor_spec, base_speed, T_peak, speed_max_ms, dt=0.1, custom_tn_df=None):
    n_max = float(motor_spec['最高轉速 (rpm)'])
    P_peak = float(motor_spec['最大功率 (kW)'])
    T_peak = float(T_peak)
    base_speed = float(base_speed)

    def get_max_torque(v):
        if v <= 0:
            v_n = 0
        else:
            v_n = v * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio
            
        if custom_tn_df is not None:
            return float(np.interp(v_n, custom_tn_df['rpm'].values, custom_tn_df['torque'].values, right=0))
            
        if v_n <= base_speed:
            return T_peak
        elif v_n <= n_max:
            return (P_peak * 1000) / (2 * math.pi * v_n / 60)
        else:
            return 0.0

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
        accels = get_safe_accel(speeds_ms, times)
        df_wltc['accel_ms2'] = accels
    
    n = len(times)
    F_roll_const = mass * G * fr
    factor_aero = 0.5 * RHO * cd * area
    wheel_circ = 2 * math.pi * wheel_radius_m
    motor_rpm = np.zeros(n)
    motor_torque = np.zeros(n)
    wheel_torque = np.zeros(n)
    
    safe_gear_eff = max(gear_eff, 0.01)
    
    for i in range(n):
        v_ms = speeds_kmh[i] / 3.6
        F_roll = F_roll_const
        F_air = factor_aero * v_ms**2
        F_accel = mass * accels[i]
        F_tractive = F_roll + F_air + F_accel
        motor_rpm[i] = v_ms / wheel_circ * 60 * gear_ratio
        
        motor_torque[i] = F_tractive * wheel_radius_m / (gear_ratio * safe_gear_eff)
        wheel_torque[i] = motor_torque[i] * gear_ratio * safe_gear_eff
        
    result_df = pd.DataFrame({
        'time': times,
        'speed_kmh': speeds_kmh,
        'accel_ms2': accels,
        'motor_rpm': motor_rpm,
        'motor_torque_Nm': motor_torque,
        'wheel_torque_Nm': wheel_torque
    })
    return result_df

# ================== 理論能耗計算函數（納入損耗對齊論文）==================
def compute_theoretical_energy_consumption(df_cycle, mass, area, cd, fr, gear_eff_percent, motor_eff_percent):
    times = df_cycle['time'].values
    speeds_kmh = df_cycle['speed_kmh'].values
    
    if 'accel_ms2' in df_cycle.columns:
        accels = df_cycle['accel_ms2'].values
    else:
        accels = get_safe_accel(speeds_kmh / 3.6, times)
    
    speeds_ms = speeds_kmh / 3.6
    
    F_roll = mass * G * fr
    F_aero = 0.5 * RHO * cd * area * speeds_ms**2
    F_acc = mass * accels
    F_trac = F_roll + F_aero + F_acc
    P_wheel = F_trac * speeds_ms
    
    eta_sys = max((gear_eff_percent / 100.0) * (motor_eff_percent / 100.0), 0.01)
    P_batt = np.zeros_like(P_wheel)
    
    drive_mask = P_wheel >= 0
    P_batt[drive_mask] = P_wheel[drive_mask] / eta_sys
    regen_mask = P_wheel < 0
    P_batt[regen_mask] = P_wheel[regen_mask] * eta_sys
    
    P_drive_batt = np.maximum(P_batt, 0)
    P_regen_batt = np.minimum(P_batt, 0)
    
    if SCIPY_AVAILABLE:
        try:
            cum_trapz_func = integrate.cumulative_trapezoid
        except AttributeError:
            cum_trapz_func = integrate.cumtrapz
            
        energy_wh_cumulative = cum_trapz_func(P_batt, times, initial=0) / 3600.0
        total_energy_wh = energy_wh_cumulative[-1]
    else:
        dt = np.diff(times, prepend=times[0])
        energy_wh_cumulative = np.cumsum(P_batt * dt) / 3600.0
        total_energy_wh = np.trapz(P_batt, times) / 3600.0

    drive_energy_wh = np.trapz(P_drive_batt, times) / 3600.0
    regen_energy_wh = np.trapz(P_regen_batt, times) / 3600.0
    wheel_energy_wh = np.trapz(P_wheel, times) / 3600.0
    distance_km = np.trapz(speeds_ms, times) / 1000.0
    
    wh_per_km_batt = float(total_energy_wh / distance_km) if distance_km > 0 else 0.0
    wh_per_km_wheel = float(wheel_energy_wh / distance_km) if distance_km > 0 else 0.0
    
    return wh_per_km_batt, wh_per_km_wheel, distance_km, total_energy_wh, drive_energy_wh, regen_energy_wh, times, P_batt, energy_wh_cumulative

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
        '放電倍率 (C)': '表示電池持續放電電流相對於容量的倍率，1C 代表可持續 1 小時放完電。',
        '串聯數': f'將多顆電池芯串聯以提高電壓。',
        '並聯數': '將多組串聯電池並聯以提高容量。',
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
st.set_page_config(layout="centered", page_title="電動載具動力估算 (WLTC 覆蓋分析)")
st.title("⚡ 電動載具動力系統估算 (WLTC 工況覆蓋分析)")

# ---------- 側邊欄（輸入參數）----------
with st.sidebar:
    st.header("🚗 整車參數規格")
    
    weight = st.number_input("車重 (kg, 不含電池)", min_value=50, value=98, step=10)
    load = st.number_input("載重 (kg)", min_value=0, value=70, step=10)
    total_mass = weight + load
    st.caption(f"總質量: {total_mass} kg")

    speed_kmh = st.number_input("目標最高車速 (km/h)", min_value=10, value=75, step=5)
    speed_ms = speed_kmh / 3.6

    st.subheader("🌬️ 阻力與環境參數設定")
    area = st.number_input("迎風面積 A (m²)", min_value=0.1, value=0.33, step=0.01, format="%.2f")
    cd = st.number_input("風阻係數 Cd", min_value=0.1, max_value=2.0, value=0.60, step=0.01, format="%.2f")
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
    avg_accel_50 = (50/3.6) / accel_time_0to50

    st.header("⛰️ 爬坡規格")
    grade_percent = st.number_input("爬坡度 (%)", min_value=0.0, value=30.0, step=0.5)

    st.header("🔧 輸入動力鍊 (Powertrain) 規格")

    with st.expander("🔹 馬達 (Motor)", expanded=True):
        voltage_option = st.radio("系統電壓", ['自動選擇', '48V', '96V'])
        voltage = None if voltage_option == '自動選擇' else int(voltage_option.replace('V', ''))
        est_mode = st.radio("估算模式", ['手動輸入', '自動估算', '讀取馬達TN曲線'], index=0)

        if est_mode == '自動估算':
            required_power, _ = calculate_power_requirements(total_mass, speed_ms, area, cd, fr)
            manual_max_power = required_power * 2
            st.info(f"⚡ 所需功率 = {required_power:.2f} kW")
            manual_peak_torque = None
            manual_max_rpm = None
            custom_tn_df = None
        elif est_mode == '手動輸入':
            manual_max_power = float(st.number_input("最大功率 (kW)", min_value=0.1, value=4.4, step=0.1))
            manual_peak_torque = float(st.number_input("最大扭矩 (Nm)", min_value=1.0, value=18.0, step=0.1))
            manual_max_rpm = float(st.number_input("最高轉速 (rpm)", min_value=100, value=9000, step=100))
            custom_tn_df = None
        else:
            tn_file = st.file_uploader("上傳 TN 曲線 (CSV)", type=["csv"], key="tn_upload")
            custom_tn_df = None
            manual_max_power, manual_peak_torque, manual_max_rpm = 4.4, 18.0, 9000.0
            if tn_file is not None:
                try:
                    df_tn = pd.read_csv(tn_file)
                    rpm_col = st.selectbox("👉 選擇轉速 (X軸) 欄位", df_tn.columns, index=0)
                    available_torque_curves = [col for col in df_tn.columns if col != rpm_col]
                    if not available_torque_curves:
                        st.error("找不到其他的扭力數據欄位！")
                    else:
                        t_col = st.selectbox("👉 選擇本次模擬使用的扭力曲線", available_torque_curves, index=0)
                        
                        df_tn_active = df_tn[[rpm_col, t_col]].copy()
                        df_tn_active[rpm_col] = pd.to_numeric(df_tn_active[rpm_col], errors='coerce')
                        df_tn_active[t_col] = pd.to_numeric(df_tn_active[t_col], errors='coerce')
                        df_tn_active = df_tn_active.dropna().sort_values(by=rpm_col)
                        
                        if len(df_tn_active) < 2:
                            st.error("有效的數值資料不足！請確保 CSV 內容為純數字。")
                        else:
                            custom_tn_df = df_tn_active.copy()
                            custom_tn_df.columns = ['rpm', 'torque']
                            custom_tn_df['power_kw'] = custom_tn_df['torque'] * custom_tn_df['rpm'] / 9550.0
                            
                            manual_max_rpm = float(custom_tn_df['rpm'].max())
                            manual_peak_torque = float(custom_tn_df['torque'].max())
                            manual_max_power = float(custom_tn_df['power_kw'].max())
                            st.info(f"📊 目前載入 [**{t_col}**]：\n最大扭矩 **{manual_peak_torque:.1f} Nm**, 最高轉速 **{manual_max_rpm:.0f} rpm**, 最大功率 **{manual_max_power:.2f} kW**")
                except Exception as e:
                    st.error(f"解析檔案失敗: {e}")

        motor_eff = st.number_input("馬達效率 (%)", min_value=1.0, max_value=100.0, value=90.0, step=1.0)

    with st.expander("🔹 齒輪 (Gear)", expanded=True):
        gear_option = st.radio("減速比", ['自動估算', '手動輸入'], index=1)
        gear_ratio = st.number_input("請輸入減速比", min_value=1.0, value=8.7, step=0.5) if gear_option == '手動輸入' else None
        gear_eff = st.number_input("齒輪箱效率 (%)", min_value=1.0, max_value=100.0, value=95.0, step=1.0)

    with st.expander("🔹 電池 (Battery)", expanded=True):
        user_battery_energy_kwh = st.number_input("電池總能量 (kWh)", min_value=0.1, value=5.0, step=0.1)
        battery_soc = st.number_input("電池可用 SOC (%)", min_value=0.0, max_value=100.0, value=90.0, step=5.0)

    # ---------- 行駛工況設定 ----------
    st.header("📁 行駛工況設定")
    wltc_file = st.file_uploader("選擇行駛工況 CSV 檔案", type=["csv"], key="wltc")
    if wltc_file is not None:
        try:
            df_wltc = pd.read_csv(wltc_file)
            time_candidates = [c for c in df_wltc.columns if 'time' in c.lower() or 't' in c.lower()]
            speed_candidates = [c for c in df_wltc.columns if 'speed' in c.lower() or 'velocity' in c.lower() or 'v' in c.lower()]
            time_col = st.selectbox("時間欄位 (秒)", df_wltc.columns, index=df_wltc.columns.get_loc(time_candidates[0]) if time_candidates else 0)
            speed_col = st.selectbox("車速欄位 (km/h)", df_wltc.columns, index=df_wltc.columns.get_loc(speed_candidates[0]) if speed_candidates else 1 if len(df_wltc.columns) > 1 else 0)
            
            df_wltc_clean = df_wltc[[time_col, speed_col]].copy()
            df_wltc_clean.columns = ['time', 'speed_kmh']
            df_wltc_clean['time'] = pd.to_numeric(df_wltc_clean['time'], errors='coerce')
            df_wltc_clean['speed_kmh'] = pd.to_numeric(df_wltc_clean['speed_kmh'], errors='coerce')
            df_wltc_clean = df_wltc_clean.dropna().sort_values(by='time')

            if len(df_wltc_clean) < 2:
                st.error("行駛工況的有效資料不足！請確保 CSV 包含數值格式的時間與車速。")
            else:
                accel_col = next((col for col in df_wltc.columns if 'accel' in col.lower() or 'a' in col.lower()), None)
                if accel_col is not None:
                    df_wltc_clean['accel_ms2'] = pd.to_numeric(df_wltc[accel_col], errors='coerce').fillna(0)
                else:
                    df_wltc_clean['accel_ms2'] = get_safe_accel(df_wltc_clean['speed_kmh'].values / 3.6, df_wltc_clean['time'].values)
                
                st.session_state.df_wltc_raw = df_wltc
                st.session_state.df_wltc_clean = df_wltc_clean
                
                gear_ratio_val = gear_ratio if gear_ratio is not None else estimate_gearbox(speed_ms, wheel_radius_m)
                df_op = compute_motor_operating_points_from_wltc(df_wltc_clean, total_mass, area, cd, fr, wheel_radius_m, gear_ratio_val, gear_eff/100.0)
                st.session_state.df_motor_operating_points = df_op
                st.success("✅ 行駛工況已載入並完成計算")
        except Exception as e:
            st.error(f"讀取檔案失敗: {e}")

# ================== 計算核心 ==================
if 'gear_ratio' not in locals() or gear_ratio is None:
    gear_ratio = estimate_gearbox(speed_ms, wheel_radius_m)
if 'gear_eff' not in locals():
    gear_eff = 95.0

if voltage is None:
    voltage = 48 if manual_max_power < 20 else 96

if est_mode == '自動估算':
    required_max_rpm = speed_ms * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio
    n_max_motor = max(required_max_rpm * 1.1, 6000.0)
else:
    n_max_motor = float(manual_max_rpm)

if est_mode == '自動估算':
    motor_spec, base_speed, T_peak = estimate_motor_from_power(manual_max_power, voltage, n_max_motor, motor_eff, base_speed=3000)
    max_power_kw_used = manual_max_power
else:
    motor_spec, base_speed, T_peak = estimate_motor_from_params(manual_max_power, manual_peak_torque, voltage, n_max_motor, motor_eff)
    max_power_kw_used = manual_max_power
    if est_mode == '讀取馬達TN曲線' and custom_tn_df is not None:
         base_speed_idx = custom_tn_df['power_kw'].argmax()
         base_speed = float(custom_tn_df.iloc[base_speed_idx]['rpm'])
         motor_spec['基速 (rpm)'] = base_speed

base_speed = float(base_speed)
T_peak = float(T_peak)

rated_power = max_power_kw_used / 2
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
    T_motor_max[over_idx], P_motor_out[over_idx] = 0, 0

# 轉換為車輪端
v_from_n = n / gear_ratio * (2 * math.pi * wheel_radius_m) * 3.6 / 60
T_wheel_max = T_motor_max * gear_ratio * ETA_DRIVE
v_max_motor = n_max_motor / gear_ratio * (2 * math.pi * wheel_radius_m) * 3.6 / 60

motor_rpm_flat, torque_flat, speed_kmh_flat, force_flat = calculate_load_curve(total_mass, area, cd, fr, wheel_radius_m, gear_ratio, speed_ms, grade_percent=0, extend_to_vmax=n_max_motor / gear_ratio * (2 * math.pi * wheel_radius_m) * 3.6 / 60)
T_wheel_flat = force_flat * wheel_radius_m

if grade_percent > 0:
    motor_rpm_climb, torque_climb, speed_kmh_climb, force_climb = calculate_load_curve(total_mass, area, cd, fr, wheel_radius_m, gear_ratio, speed_ms, grade_percent, extend_to_vmax=n_max_motor / gear_ratio * (2 * math.pi * wheel_radius_m) * 3.6 / 60)
    T_wheel_climb = force_climb * wheel_radius_m
else:
    motor_rpm_climb, torque_climb = None, None

time_acc, speed_acc, disp_acc = simulate_acceleration(total_mass, area, cd, fr, wheel_radius_m, gear_ratio, motor_spec, base_speed, T_peak, speed_ms, dt=0.1, custom_tn_df=custom_tn_df)
actual_0to50 = time_acc[np.argmax(speed_acc >= 50)] if np.any(speed_acc >= 50) else np.inf
actual_full_time = time_acc[np.argmax(speed_acc >= speed_kmh * 0.99)] if np.any(speed_acc >= speed_kmh * 0.99) else np.inf

if "default_motor_spec" not in st.session_state: st.session_state.default_motor_spec = motor_spec.copy()
if "default_battery_spec" not in st.session_state: st.session_state.default_battery_spec = battery_spec.copy()

# ================== 顯示區 ==================
st.subheader("📋 規格摘要")
with st.expander("📦 馬達規格", expanded=True):
    st.markdown(render_json_with_diff(motor_spec, st.session_state.default_motor_spec), unsafe_allow_html=True)
with st.expander("🏎️ 動態性能表現", expanded=True):
    st.metric("0→最高車速起步所需馬達扭矩", f"{T_motor_start_full:.1f} Nm")
    st.metric("0→50 km/h 起步所需馬達扭矩", f"{T_motor_start_50:.1f} Nm")
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

# ================== 圖1：馬達 TN 曲線 + 功率曲線 + 工作點 ==================
st.markdown("## 📈 圖1：馬達 TN 曲線 + 功率曲線 + 工作點")
x_upper = float(n_max_motor * 1.1)
grid_step = float(T_peak / 4.0 if T_peak > 0 else 10)
y_min_raw = min(0, float(T_motor_max.min()), float(torque_flat.min()))
if "df_motor_operating_points" in st.session_state:
    min_op = st.session_state.df_motor_operating_points['motor_torque_Nm'].min()
    if not np.isnan(min_op): y_min_raw = min(y_min_raw, float(min_op))
y_min_torque = math.floor(y_min_raw / grid_step) * grid_step

y_max_torque = T_peak + grid_step
if "df_motor_operating_points" in st.session_state:
    max_op = st.session_state.df_motor_operating_points['motor_torque_Nm'].max()
    if not np.isnan(max_op) and max_op > y_max_torque: y_max_torque = math.ceil(max_op / grid_step) * grid_step

ratio = float(max_power_kw_used / T_peak if T_peak > 0 else 1)
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
    fig1.add_trace(go.Scatter(x=motor_rpm_climb, y=torque_climb, mode='lines', name='爬坡負載線', line=dict(color='green', width=3, dash='dot')), secondary_y=False)
fig1.add_trace(go.Scatter(x=n, y=P_motor_out, mode='lines', name='馬達功率', line=dict(color='gold', width=2, dash='solid')), secondary_y=True)
if "df_motor_operating_points" in st.session_state:
    fig1.add_trace(go.Scatter(x=st.session_state.df_motor_operating_points['motor_rpm'], y=st.session_state.df_motor_operating_points['motor_torque_Nm'], mode='markers', marker=dict(size=4, color='cyan', opacity=0.6), name='工作點'), secondary_y=False)

fig1.update_yaxes(title_text="扭矩 (Nm)", secondary_y=False, range=[y_min_torque, y_max_torque], tickvals=y_ticks, tickfont=dict(color='white'), zeroline=True, zerolinecolor='gray')
fig1.update_yaxes(title_text="功率 (kW)", secondary_y=True, range=[p_min, p_max], tickvals=p_ticks, tickfont=dict(color='white'), showgrid=False)
fig1.update_xaxes(title_text="轉速 (rpm)", range=[0, x_upper], tickvals=x_ticks, tickfont=dict(color='white'))
fig1.add_annotation(x=base_speed, y=T_peak, text=f"基速: {base_speed:.0f}", showarrow=True, arrowhead=2, arrowcolor="green", ax=0, ay=-45, font=dict(color="lightgreen", size=12), bgcolor="rgba(26,28,35,0.9)", bordercolor="green", borderwidth=1)

# ================= 新增：交點標記與防重疊數值標籤 (圖1) =================
intersections_flat_motor = find_intersection(n, T_motor_max, motor_rpm_flat, torque_flat)
for i, (x_cross, y_cross) in enumerate(intersections_flat_motor):
    fig1.add_trace(go.Scatter(x=[x_cross], y=[y_cross], mode='markers', name='平路交點' if i==0 else None, marker=dict(color='red', size=12, symbol='x'), showlegend=(i==0)), secondary_y=False)
    fig1.add_annotation(x=x_cross, y=y_cross, text=f'{x_cross:.0f} rpm<br>{y_cross:.1f} Nm', showarrow=True, arrowhead=2, ax=45, ay=40, font=dict(size=11, color="white"), bgcolor="rgba(255,0,0,0.4)", bordercolor="red", borderwidth=1)

if motor_rpm_climb is not None:
    intersections_climb_motor = find_intersection(n, T_motor_max, motor_rpm_climb, torque_climb)
    for i, (x_cross, y_cross) in enumerate(intersections_climb_motor):
        fig1.add_trace(go.Scatter(x=[x_cross], y=[y_cross], mode='markers', name='爬坡交點' if i==0 else None, marker=dict(color='green', size=12, symbol='x'), showlegend=(i==0)), secondary_y=False)
        fig1.add_annotation(x=x_cross, y=y_cross, text=f'{x_cross:.0f} rpm<br>{y_cross:.1f} Nm', showarrow=True, arrowhead=2, ax=-45, ay=-40, font=dict(size=11, color="white"), bgcolor="rgba(0,128,0,0.4)", bordercolor="green", borderwidth=1)

fig1.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5), margin=dict(l=80, r=110, t=100, b=20), height=550)
st.plotly_chart(fig1, use_container_width=True)

st.markdown("---")
st.markdown("## 📈 圖2：車輪扭矩 vs 車速 + 工作點")
fig2 = go.Figure()
fig2.add_trace(go.Scatter(x=v_from_n, y=T_wheel_max, mode='lines', name='最大車輪扭矩', line=dict(color='dodgerblue', width=3)))
fig2.add_trace(go.Scatter(x=speed_kmh_flat, y=T_wheel_flat, mode='lines', name='平路負載線', line=dict(color='red', width=3, dash='dash')))
if T_wheel_climb is not None:
    fig2.add_trace(go.Scatter(x=speed_kmh_climb, y=T_wheel_climb, mode='lines', name='爬坡負載線', line=dict(color='green', width=3, dash='dot')))
if "df_motor_operating_points" in st.session_state:
    fig2.add_trace(go.Scatter(x=st.session_state.df_motor_operating_points['speed_kmh'], y=st.session_state.df_motor_operating_points['wheel_torque_Nm'], mode='markers', marker=dict(size=4, color='cyan', opacity=0.6), name='工作點'))

T_design_flat = float(T_wheel_flat[np.argmin(np.abs(speed_kmh_flat - speed_kmh))])
T_at_vmax = float(np.interp(v_max_motor, v_from_n, T_wheel_max) if v_max_motor <= v_from_n.max() else 0)
T_wheel_peak = float(T_wheel_max.max())

y_min_raw_w = min(0, float(T_wheel_max.min()), float(T_wheel_flat.min()))
if "df_motor_operating_points" in st.session_state:
    min_op_w = st.session_state.df_motor_operating_points['wheel_torque_Nm'].min()
    if not np.isnan(min_op_w): y_min_raw_w = min(y_min_raw_w, float(min_op_w))

grid_step_wheel = T_wheel_peak / 4.0 if T_wheel_peak > 0 else 10
y_min_wheel = math.floor(y_min_raw_w / grid_step_wheel) * grid_step_wheel
y_max_wheel = T_wheel_peak + grid_step_wheel

fig2.update_yaxes(title_text="車輪扭矩 (Nm)", range=[y_min_wheel, y_max_wheel], tickfont=dict(color='white'), zeroline=True, zerolinecolor='gray')
fig2.update_xaxes(title_text="車速 (km/h)", range=[0, max(v_max_motor, speed_kmh) * 1.15 if max(v_max_motor, speed_kmh) > 0 else 100], tickfont=dict(color='white'), zeroline=True, zerolinecolor='gray')
fig2.add_annotation(x=speed_kmh, y=T_design_flat, text=f"<b>目標: {speed_kmh:.0f} km/h</b>", showarrow=True, arrowhead=2, arrowcolor="orange", ax=-55, ay=-65, font=dict(color="orange", size=12), bgcolor="rgba(26,28,35,0.9)", bordercolor="orange", borderwidth=1)
fig2.add_annotation(x=v_max_motor, y=T_at_vmax, text=f"<b>極速: {v_max_motor:.0f} km/h</b>", showarrow=True, arrowhead=2, arrowcolor="purple", ax=65, ay=-45, font=dict(color="#d8b4e2", size=12), bgcolor="rgba(26,28,35,0.9)", bordercolor="purple", borderwidth=1)

# ================= 新增：交點標記與防重疊數值標籤 (圖2) =================
intersections_flat_wheel = find_intersection(v_from_n, T_wheel_max, speed_kmh_flat, T_wheel_flat)
for i, (x_cross, y_cross) in enumerate(intersections_flat_wheel):
    fig2.add_trace(go.Scatter(x=[x_cross], y=[y_cross], mode='markers', name='平路交點' if i==0 else None, marker=dict(color='red', size=12, symbol='x'), showlegend=(i==0)))
    fig2.add_annotation(x=x_cross, y=y_cross, text=f'{x_cross:.1f} km/h<br>{y_cross:.1f} Nm', showarrow=True, arrowhead=2, ax=55, ay=50, font=dict(size=11, color="white"), bgcolor="rgba(255,0,0,0.4)", bordercolor="red", borderwidth=1)

if T_wheel_climb is not None:
    intersections_climb_wheel = find_intersection(v_from_n, T_wheel_max, speed_kmh_climb, T_wheel_climb)
    for i, (x_cross, y_cross) in enumerate(intersections_climb_wheel):
        fig2.add_trace(go.Scatter(x=[x_cross], y=[y_cross], mode='markers', name='爬坡交點' if i==0 else None, marker=dict(color='green', size=12, symbol='x'), showlegend=(i==0)))
        fig2.add_annotation(x=x_cross, y=y_cross, text=f'{x_cross:.1f} km/h<br>{y_cross:.1f} Nm', showarrow=True, arrowhead=2, ax=-55, ay=50, font=dict(size=11, color="white"), bgcolor="rgba(0,128,0,0.4)", bordercolor="green", borderwidth=1)


fig2.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5), margin=dict(l=80, r=110, t=100, b=20), height=550)
st.plotly_chart(fig2, use_container_width=True)

st.markdown("---")
st.markdown("## 📈 圖3：加速性能（速度與位移 vs 時間）")
fig3 = make_subplots(specs=[[{"secondary_y": True}]])
fig3.add_trace(go.Scatter(x=time_acc, y=speed_acc, mode='lines', name='車速', line=dict(color='dodgerblue', width=3)), secondary_y=False)
fig3.add_trace(go.Scatter(x=time_acc, y=disp_acc, mode='lines', name='位移', line=dict(color='red', width=2, dash='dash')), secondary_y=True)
fig3.update_layout(height=400, margin=dict(l=20, r=20, t=40, b=20))
st.plotly_chart(fig3, use_container_width=True)
st.markdown("---")

if "df_wltc_clean" in st.session_state:
    st.markdown("## 📈 圖4：行駛工況分析")
    df_p = st.session_state.df_wltc_clean
    fig4 = make_subplots(specs=[[{"secondary_y": True}]])
    fig4.add_trace(go.Scatter(x=df_p['time'], y=df_p['speed_kmh'], name='車速', line=dict(color='dodgerblue')), secondary_y=False)
    fig4.add_trace(go.Scatter(x=df_p['time'], y=df_p['accel_ms2'], name='加速度', line=dict(color='red', dash='dash')), secondary_y=True)
    fig4.update_layout(height=400, margin=dict(l=20, r=20, t=40, b=20))
    st.plotly_chart(fig4, use_container_width=True)
    st.markdown("---")

# ================== 圖5：理論能耗分析 & 里程估算 ==================
if "df_wltc_clean" in st.session_state:
    st.markdown("## 📈 圖5：理論能耗分析與里程預估 (含損耗對齊)")
    st.caption("本模組基於行駛工況積分計算能耗，並根據側邊欄設定的**電池容量與 SOC** 進行續航里程預估。")
    
    df_energy = st.session_state.df_wltc_clean.copy()
    gear_eff_val = gear_eff if 'gear_eff' in locals() else 95.0
    motor_eff_val = motor_eff if 'motor_eff' in locals() else 90.0
    
    wh_per_km_batt, wh_per_km_wheel, total_dist_km, total_energy_wh, drive_energy_wh, regen_energy_wh, times, power_batt_w, energy_cum = compute_theoretical_energy_consumption(
        df_energy, total_mass, area, cd, fr, gear_eff_val, motor_eff_val
    )
    
    usable_energy_wh = user_battery_energy_kwh * 1000 * (battery_soc / 100.0)
    estimated_range_km = usable_energy_wh / wh_per_km_batt if wh_per_km_batt > 0 else 0
    
    c1, c2, c3 = st.columns(3)
    c1.metric("🏁 工況總里程", f"{total_dist_km:.3f} km")
    c2.metric("⚙️ 輪上淨能耗 (純物理)", f"{wh_per_km_wheel:.2f} Wh/km")
    c3.metric("🔋 電池理論能耗 (含損耗)", f"{wh_per_km_batt:.2f} Wh/km")
    
    st.markdown("#### 🔋 續航里程估算結果")
    cc1, cc2, cc3 = st.columns(3)
    cc1.metric("⚡ 可用總能量", f"{usable_energy_wh:.1f} Wh")
    cc2.metric("🎯 理論預估里程", f"{estimated_range_km:.1f} km", delta=f"基於 {user_battery_energy_kwh:.1f}kWh 電池")
    cc3.metric("📈 驅動耗電 (電池端)", f"{drive_energy_wh:.1f} Wh")
    
    fig5 = make_subplots(specs=[[{"secondary_y": True}]])
    fig5.add_trace(go.Scatter(x=times, y=np.maximum(power_batt_w, 0), fill='tozeroy', name='驅動輸出', line=dict(color='crimson', width=1)), secondary_y=False)
    fig5.add_trace(go.Scatter(x=times, y=np.minimum(power_batt_w, 0), fill='tozeroy', name='回收充電', line=dict(color='seagreen', width=1)), secondary_y=False)
    fig5.add_trace(go.Scatter(x=times, y=energy_cum, name='累積淨耗電', line=dict(color='gold', width=3)), secondary_y=True)
    
    fig5.update_layout(height=450, margin=dict(l=20, r=20, t=40, b=20), legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5))
    st.plotly_chart(fig5, use_container_width=True)