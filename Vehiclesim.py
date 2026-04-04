import streamlit as st
import pandas as pd
import numpy as np
import math
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from io import BytesIO

# ================== Constants & Default Parameters ==================
G = 9.81
RHO = 1.2
FR = 0.015
CD = 0.4
ETA_DRIVE = 0.9                # Drivetrain efficiency
ETA_MOTOR = 1.0                 # Motor efficiency (set to 1, not affecting calculations)
ETA_CONTROLLER = 0.95
BATTERY_ENERGY_DENSITY = 150    # Wh/kg
MOTOR_POWER_DENSITY = 1.0       # kW/kg
CELL_VOLTAGE = 3.7              # Li-ion cell voltage (V)
CELL_CAPACITY = 2.5             # Ah per cell (18650 typical)

# ================== Helper Functions ==================
def get_cd_by_vehicle(vehicle_type):
    mapping = {
        'Small EV': 0.3,
        'Electric Scooter': 0.5,
        'Electric Tricycle': 0.45,
        'Golf Cart': 0.4
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
        'Type': 'PMSM (Permanent Magnet Synchronous Motor)',
        'Max Power (kW)': round(max_power_kw, 2),
        'Rated Power (kW)': round(rated_power, 2),
        'Peak Torque (Nm)': round(T_peak, 1),
        'Rated Torque (Nm)': round(T_rated, 1),
        'Rated Speed (rpm)': base_speed,
        'Max Speed (rpm)': round(n_max, 0),
        'Base Speed (rpm)': base_speed,
        'Voltage (V)': voltage,
        'Efficiency (%)': round(ETA_MOTOR * 100, 1),
        'Est. Weight (kg)': round(max_power_kw * MOTOR_POWER_DENSITY, 1)
    }
    return motor_spec, base_speed, T_peak

def estimate_motor_from_params(max_power_kw, peak_torque_Nm, voltage, n_max):
    base_speed = (max_power_kw * 1000 * 60) / (2 * math.pi * peak_torque_Nm)
    if base_speed > n_max:
        base_speed = n_max
    rated_power = max_power_kw / 2
    T_rated = peak_torque_Nm / 2
    motor_spec = {
        'Type': 'PMSM (Permanent Magnet Synchronous Motor)',
        'Max Power (kW)': round(max_power_kw, 2),
        'Rated Power (kW)': round(rated_power, 2),
        'Peak Torque (Nm)': round(peak_torque_Nm, 1),
        'Rated Torque (Nm)': round(T_rated, 1),
        'Rated Speed (rpm)': round(base_speed, 0),
        'Max Speed (rpm)': round(n_max, 0),
        'Base Speed (rpm)': round(base_speed, 0),
        'Voltage (V)': voltage,
        'Efficiency (%)': round(ETA_MOTOR * 100, 1),
        'Est. Weight (kg)': round(max_power_kw * MOTOR_POWER_DENSITY, 1)
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
        'Type': 'Li-ion Battery',
        'Nominal Voltage (V)': voltage,
        'Capacity (Ah)': round(capacity_ah, 1),
        'Energy (kWh)': round(energy_kwh, 2),
        'Discharge Rate (C)': round(c_rate, 1),
        'Series Cells': series,
        'Parallel Cells': parallel,
        'Est. Weight (kg)': round(weight, 1)
    }
    return battery_spec

def estimate_controller(max_power_kw, voltage):
    I_max = (max_power_kw * 1000) / voltage
    controller_spec = {
        'Type': 'MOSFET Controller',
        'Max Power (kW)': round(max_power_kw, 2),
        'Voltage Range (V)': f"{int(voltage*0.8)}-{int(voltage*1.2)}",
        'Max Current (A)': round(I_max, 1),
        'Efficiency (%)': ETA_CONTROLLER * 100
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
    n_max = motor_spec['Max Speed (rpm)']
    P_peak = motor_spec['Max Power (kW)']

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

# ================== Streamlit Interface ==================
st.set_page_config(layout="centered", page_title="EV Powertrain Estimator (Mobile Optimized)")

st.title("⚡ EV Powertrain Estimator (Mobile Optimized)")

# ---------- Sidebar (Input Parameters) ----------
with st.sidebar:
    st.header("🚗 Vehicle Parameters")
    vehicle_type = st.selectbox("Vehicle Type", ['Small EV', 'Electric Scooter', 'Electric Tricycle', 'Golf Cart'], index=1)
    weight = st.number_input("Vehicle Weight (kg, excluding battery)", min_value=50, value=98, step=10)
    load = st.number_input("Payload (kg)", min_value=0, value=63, step=10)
    total_mass = weight + load
    st.caption(f"Gross Mass: {total_mass} kg")

    speed_kmh = st.number_input("Target Top Speed (km/h)", min_value=10, value=75, step=5)
    speed_ms = speed_kmh / 3.6

    area = st.number_input("Frontal Area (m²)", min_value=0.3, value=0.61, step=0.05, format="%.2f")

    st.subheader("Tire Specification")
    tire_width = st.number_input("Tire Width (mm)", min_value=50, value=110, step=5)
    tire_aspect = st.number_input("Aspect Ratio (%)", min_value=30, value=70, step=5)
    rim_dia_inch = st.number_input("Rim Diameter (inch)", min_value=8, value=12, step=1,
                                   help="This is the rim diameter (wheel inner diameter) used to calculate tire radius.")
    st.caption("Note: This is the rim diameter (inner diameter), not the outer tire radius.")

    sidewall_height_mm = tire_width * tire_aspect / 100
    rim_radius_mm = (rim_dia_inch * 25.4) / 2
    tire_radius_m = (rim_radius_mm + sidewall_height_mm) / 1000
    st.caption(f"Calculated Tire Radius: {tire_radius_m:.4f} m")
    wheel_radius_m = tire_radius_m

    voltage_option = st.radio("System Voltage", ['Auto', '48V', '96V'])
    if voltage_option == 'Auto':
        voltage = None
    else:
        voltage = int(voltage_option.replace('V', ''))

    gear_option = st.radio("Gear Ratio", ['Auto Estimate', 'Manual Input'])
    if gear_option == 'Manual Input':
        gear_ratio = st.number_input("Enter Gear Ratio", min_value=1.0, value=8.7, step=0.5)
    else:
        gear_ratio = None

    # ---------- Motor Sizing ----------
    st.markdown("---")
    st.subheader("🔧 Motor Sizing")
    est_mode = st.radio("Estimation Mode", ['Auto', 'Manual'], index=0,
                        help="Auto: Calculates required power based on target speed. Manual: You specify max power and torque independently.")

    if est_mode == 'Auto':
        cd = get_cd_by_vehicle(vehicle_type)
        fr = FR
        required_power, _ = calculate_power_requirements(total_mass, speed_ms, area, cd, fr)
        max_power_kw = required_power * 2
        st.info(f"⚡ Required power = {required_power:.2f} kW → Max Power = {max_power_kw:.2f} kW")
        manual_max_power = max_power_kw
        manual_peak_torque = None
    else:
        manual_max_power = st.number_input("Max Power (kW)", min_value=0.1, value=10.24, step=0.1)
        manual_peak_torque = st.number_input("Peak Torque (Nm)", min_value=1.0, value=32.6, step=0.1)
        base_speed_calc = (manual_max_power * 1000 * 60) / (2 * math.pi * manual_peak_torque)
        st.caption(f"Corresponding Base Speed ≈ {base_speed_calc:.0f} rpm")

    # ---------- Acceleration Targets ----------
    st.markdown("---")
    st.subheader("⚡ Acceleration Targets")
    accel_time_full = st.number_input("0→Top Speed Time (s)", min_value=1.0, value=10.0, step=0.5)
    avg_accel_full = speed_ms / accel_time_full

    accel_time_0to50 = st.number_input("0→50 km/h Time (s)", min_value=1.0, value=5.0, step=0.5)
    speed_50_ms = 50 / 3.6
    avg_accel_50 = speed_50_ms / accel_time_0to50

    # ---------- Gradeability ----------
    st.markdown("---")
    st.subheader("⛰️ Gradeability")
    grade_percent = st.number_input("Grade (%)", min_value=0.0, value=0.0, step=0.5)

    use_range = st.checkbox("Specify Range")
    if use_range:
        desired_range = st.number_input("Desired Range (km)", min_value=1, value=50, step=5)
    else:
        desired_range = None

    st.markdown("---")
    st.caption("Modify parameters and results update automatically below.")

# ================== Core Calculations ==================
cd = get_cd_by_vehicle(vehicle_type)
fr = FR

if voltage is None:
    if est_mode == 'Auto':
        power_val = manual_max_power
    else:
        power_val = manual_max_power
    voltage = 48 if power_val < 20 else 96

if gear_ratio is None:
    gear_ratio = estimate_gearbox(speed_ms, wheel_radius_m)

required_max_rpm = speed_ms * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio
n_max_motor = max(required_max_rpm * 1.1, 6000)

if est_mode == 'Auto':
    motor_spec, base_speed, T_peak = estimate_motor_from_power(manual_max_power, voltage, n_max_motor, base_speed=3000)
    max_power_kw_used = manual_max_power
else:
    motor_spec, base_speed, T_peak = estimate_motor_from_params(manual_max_power, manual_peak_torque, voltage, n_max_motor)
    max_power_kw_used = manual_max_power

rated_power = max_power_kw_used / 2

if desired_range:
    avg_speed_ms = speed_ms * 0.7
    avg_speed_kmh = avg_speed_ms * 3.6
    time_h = desired_range / avg_speed_kmh
    avg_power_kw = rated_power * 0.7
    battery_spec = estimate_battery(avg_power_kw, voltage, duration_h=time_h)
else:
    battery_spec = estimate_battery(rated_power, voltage, duration_h=1.0)

controller_spec = estimate_controller(max_power_kw_used, voltage)

gearbox_spec = {
    'Type': 'Fixed Ratio Gearbox',
    'Gear Ratio': round(gear_ratio, 2),
    'Efficiency (%)': 95
}

# ---------- Launch Torque Requirements ----------
F_roll_start = total_mass * G * fr
F_accel_full = total_mass * avg_accel_full
F_total_start_full = F_roll_start + F_accel_full
T_wheel_start_full = F_total_start_full * wheel_radius_m
T_motor_start_full = T_wheel_start_full / (gear_ratio * ETA_DRIVE)

F_accel_50 = total_mass * avg_accel_50
F_total_start_50 = F_roll_start + F_accel_50
T_wheel_start_50 = F_total_start_50 * wheel_radius_m
T_motor_start_50 = T_wheel_start_50 / (gear_ratio * ETA_DRIVE)

# ---------- Load Lines ----------
motor_rpm_flat, torque_flat, speed_kmh_flat, force_flat = calculate_load_curve(
    total_mass, area, cd, fr, wheel_radius_m, gear_ratio, speed_ms, grade_percent=0
)

if grade_percent > 0:
    motor_rpm_climb, torque_climb, speed_kmh_climb, force_climb = calculate_load_curve(
        total_mass, area, cd, fr, wheel_radius_m, gear_ratio, speed_ms, grade_percent
    )
else:
    motor_rpm_climb, torque_climb = None, None

# ---------- 0-50 Acceleration Demand Curve ----------
F_accel_const = total_mass * avg_accel_50
T_accel_const_motor = F_accel_const * wheel_radius_m / (gear_ratio * ETA_DRIVE)
v_accel = np.linspace(0, min(50, speed_kmh), 50)
torque_flat_at_v = np.interp(v_accel, speed_kmh_flat, torque_flat)
torque_total_accel_motor = torque_flat_at_v + T_accel_const_motor
torque_total_accel_wheel = torque_total_accel_motor * gear_ratio * ETA_DRIVE

# ---------- Acceleration Simulation ----------
time_acc, speed_acc, disp_acc = simulate_acceleration(
    total_mass, area, cd, fr, wheel_radius_m, gear_ratio,
    motor_spec, base_speed, T_peak, speed_ms, dt=0.1
)

# Actual 0-50 time
if np.any(speed_acc >= 50):
    actual_0to50 = time_acc[np.argmax(speed_acc >= 50)]
else:
    actual_0to50 = np.inf

# Actual 0-top speed time
if np.any(speed_acc >= speed_kmh * 0.99):
    actual_full_time = time_acc[np.argmax(speed_acc >= speed_kmh * 0.99)]
else:
    actual_full_time = np.inf

# ================== Display Area (Single Column Vertical) ==================

# ---------- Summary ----------
st.subheader("📋 Specification Summary")

with st.expander("📦 Motor Specifications", expanded=True):
    st.json(motor_spec)

with st.expander("🚦 Launch Performance Comparison", expanded=True):
    st.metric("Required Motor Torque (0→Top Speed)", f"{T_motor_start_full:.1f} Nm")
    st.metric("Required Motor Torque (0→50 km/h)", f"{T_motor_start_50:.1f} Nm")
    if T_motor_start_full <= T_peak and T_motor_start_50 <= T_peak:
        st.success("✅ Motor peak torque sufficient for both acceleration demands")
    else:
        short = max(0, T_motor_start_full - T_peak, T_motor_start_50 - T_peak)
        st.error(f"❌ Motor peak torque insufficient, need +{short:.1f} Nm")

with st.expander("⚡ Acceleration Performance Comparison", expanded=True):
    col_a, col_b = st.columns(2)
    with col_a:
        st.metric("Target 0→50 km/h", f"{accel_time_0to50:.1f} s")
        st.metric("Actual 0→50 km/h", f"{actual_0to50:.1f} s")
        if actual_0to50 <= accel_time_0to50:
            st.success("✅ Target met")
        else:
            st.error("❌ Not met")
    with col_b:
        st.metric("Target 0→Top Speed", f"{accel_time_full:.1f} s")
        st.metric("Actual 0→Top Speed", f"{actual_full_time:.1f} s")
        if actual_full_time <= accel_time_full:
            st.success("✅ Target met")
        else:
            st.error("❌ Not met")

with st.expander("🔋 Battery", expanded=False):
    st.json(battery_spec)

with st.expander("🎛️ Controller", expanded=False):
    st.json(controller_spec)

with st.expander("⚙️ Gearbox", expanded=False):
    st.json(gearbox_spec)

with st.expander("🔁 Conversion Factors", expanded=False):
    torque_factor = gear_ratio * ETA_DRIVE
    speed_factor = (2 * math.pi * wheel_radius_m * 60) / (gear_ratio * 1000) * 3.6
    st.metric("Wheel Torque / Motor Torque", f"{torque_factor:.3f}")
    st.caption("Formula: gear ratio × drivetrain efficiency")
    st.metric("Speed (km/h) / Motor Speed (rpm)", f"{speed_factor:.6f}")
    st.caption("Formula: (2π × tire radius(m) × 60) / (gear ratio × 1000) × 3.6")

idx_design_local = np.argmin(np.abs(speed_kmh_flat - speed_kmh))
T_design_flat_local = (force_flat[idx_design_local] * wheel_radius_m)
F_design_flat_local = force_flat[idx_design_local]
with st.expander("🔧 Performance at Design Top Speed", expanded=False):
    st.metric("Wheel Torque at Top Speed", f"{T_design_flat_local:.1f} Nm")
    st.metric("Wheel Thrust at Top Speed", f"{F_design_flat_local:.1f} N")

# Download Excel
df_motor = pd.DataFrame([motor_spec])
df_battery = pd.DataFrame([battery_spec])
df_controller = pd.DataFrame([controller_spec])
df_gearbox = pd.DataFrame([gearbox_spec])

output = BytesIO()
with pd.ExcelWriter(output, engine='openpyxl') as writer:
    df_motor.to_excel(writer, sheet_name='Motor', index=False)
    df_battery.to_excel(writer, sheet_name='Battery', index=False)
    df_controller.to_excel(writer, sheet_name='Controller', index=False)
    df_gearbox.to_excel(writer, sheet_name='Gearbox', index=False)
st.download_button(
    label="📥 Download Excel Report",
    data=output.getvalue(),
    file_name="powertrain_spec.xlsx",
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    use_container_width=True
)

st.markdown("---")

# ---------- Figure 1: Motor TN Curve + Power Curve ----------
n_max_motor = motor_spec['Max Speed (rpm)']
P_peak = motor_spec['Max Power (kW)']

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
    go.Scatter(x=n, y=T_motor_max, mode='lines', name='Motor Max Torque', line=dict(color='blue', width=3)),
    secondary_y=False
)
fig1.add_trace(
    go.Scatter(x=motor_rpm_flat, y=torque_flat, mode='lines', name='Flat Road Load (Motor Side)', line=dict(color='red', width=3, dash='dash')),
    secondary_y=False
)
if motor_rpm_climb is not None:
    fig1.add_trace(
        go.Scatter(x=motor_rpm_climb, y=torque_climb, mode='lines', name=f'Grade Load ({grade_percent}%)',
                   line=dict(color='green', width=3, dash='dot')),
        secondary_y=False
    )
fig1.add_trace(
    go.Scatter(x=n, y=P_motor_out, mode='lines', name='Motor Power', line=dict(color='gold', width=2, dash='solid')),
    secondary_y=True
)

# Key points
fig1.add_trace(
    go.Scatter(x=[0], y=[T_peak], mode='markers+text', name='Peak Torque Point',
               text=[f'{T_peak:.1f} Nm'], textposition='bottom right',
               marker=dict(color='blue', size=10), textfont=dict(size=10)),
    secondary_y=False
)
fig1.add_trace(
    go.Scatter(x=[base_speed], y=[T_peak], mode='markers+text', name='Base Speed Point',
               text=[f'Base Speed: {base_speed:.0f} rpm'], textposition='top left',
               marker=dict(color='green', size=10), textfont=dict(size=10)),
    secondary_y=False
)
T_at_max_n = (P_peak * 1000) / (2 * math.pi * n_max_motor / 60) if n_max_motor > 0 else 0
fig1.add_trace(
    go.Scatter(x=[n_max_motor], y=[T_at_max_n], mode='markers+text', name='Max Speed Point',
               text=[f'{n_max_motor:.0f} rpm, {T_at_max_n:.1f} Nm'],
               textposition='top right',
               marker=dict(color='purple', size=10), textfont=dict(size=10)),
    secondary_y=False
)
fig1.add_vline(x=design_rpm, line_width=2, line_dash="dash", line_color="orange", opacity=0.9)
T_at_design = np.interp(design_rpm, n, T_motor_max) if design_rpm <= n_max_motor else 0
fig1.add_trace(
    go.Scatter(x=[design_rpm], y=[T_at_design], mode='markers+text',
               name='Design Speed Corresponding RPM',
               text=[f'{design_rpm:.0f} rpm, {T_at_design:.1f} Nm'],
               textposition='top center',
               marker=dict(color='orange', size=10),
               textfont=dict(size=10)),
    secondary_y=False
)

# Intersections
intersections_flat = find_intersection(n, T_motor_max, motor_rpm_flat, torque_flat)
for i, (x_cross, y_cross) in enumerate(intersections_flat):
    fig1.add_trace(
        go.Scatter(x=[x_cross], y=[y_cross], mode='markers',
                   name=f'Flat Road Intersection {i+1}' if i==0 else None,
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
                       name=f'Grade Intersection {i+1}' if i==0 else None,
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
fig1.update_xaxes(title_text="Speed (rpm)")
fig1.update_yaxes(title_text="Torque (Nm)", secondary_y=False)
fig1.update_yaxes(title_text="Power (kW)", secondary_y=True)
st.plotly_chart(fig1, use_container_width=True)

st.markdown("---")

# ---------- Figure 2: Wheel Torque vs Vehicle Speed ----------
v_from_n = n / gear_ratio * (2 * math.pi * wheel_radius_m) * 3.6 / 60
T_wheel_max = T_motor_max * gear_ratio * ETA_DRIVE
T_wheel_flat = force_flat * wheel_radius_m

if grade_percent > 0:
    T_wheel_climb = force_climb * wheel_radius_m
    v_climb = speed_kmh_climb
else:
    T_wheel_climb = None
    v_climb = None

idx_design = np.argmin(np.abs(speed_kmh_flat - speed_kmh))
T_design_flat = T_wheel_flat[idx_design]
v_max_motor = n_max_motor / gear_ratio * (2 * math.pi * wheel_radius_m) * 3.6 / 60
T_at_vmax = np.interp(v_max_motor, v_from_n, T_wheel_max) if v_max_motor <= v_from_n.max() else 0

fig2 = go.Figure()
fig2.add_trace(go.Scatter(x=v_from_n, y=T_wheel_max, mode='lines', name='Max Wheel Torque',
                           line=dict(color='blue', width=3)))
fig2.add_trace(go.Scatter(x=speed_kmh_flat, y=T_wheel_flat, mode='lines', name='Flat Road Load (Wheel)',
                           line=dict(color='red', width=3, dash='dash')))
if T_wheel_climb is not None:
    fig2.add_trace(go.Scatter(x=v_climb, y=T_wheel_climb, mode='lines',
                               name=f'Grade Load ({grade_percent}%)',
                               line=dict(color='green', width=3, dash='dot')))
fig2.add_trace(go.Scatter(x=v_accel, y=torque_total_accel_wheel, mode='lines',
                           name=f'0-50km/h Acceleration Demand ({accel_time_0to50}s)',
                           line=dict(color='blue', width=2, dash='dash')))

fig2.add_vline(x=speed_kmh, line_width=2, line_dash="dash", line_color="orange", opacity=0.9)
fig2.add_trace(go.Scatter(x=[speed_kmh], y=[T_design_flat], mode='markers+text',
                           name='Design Top Speed Point',
                           text=[f'{speed_kmh:.0f} km/h, {T_design_flat:.1f} Nm'],
                           textposition='top right',
                           marker=dict(color='orange', size=10),
                           textfont=dict(size=10)))
fig2.add_trace(go.Scatter(x=[v_max_motor], y=[T_at_vmax], mode='markers+text',
                           name='Motor Max Speed Corresponding Speed',
                           text=[f'Max Motor Speed\n{v_max_motor:.0f} km/h, {T_at_vmax:.1f} Nm'],
                           textposition='top left',
                           marker=dict(color='purple', size=10),
                           textfont=dict(size=10)))

intersections_flat_wheel = find_intersection(v_from_n, T_wheel_max, speed_kmh_flat, T_wheel_flat)
for i, (x_cross, y_cross) in enumerate(intersections_flat_wheel):
    fig2.add_trace(
        go.Scatter(x=[x_cross], y=[y_cross], mode='markers',
                   name=f'Flat Road Intersection {i+1}' if i==0 else None,
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
                       name=f'Grade Intersection {i+1}' if i==0 else None,
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
fig2.update_xaxes(title_text="Vehicle Speed (km/h)", range=[0, x_max])
fig2.update_yaxes(title_text="Torque (Nm)")
st.plotly_chart(fig2, use_container_width=True)

st.markdown("---")

# ---------- Figure 3: Speed & Displacement vs Time ----------
fig3 = make_subplots(specs=[[{"secondary_y": True}]])
fig3.add_trace(
    go.Scatter(x=time_acc, y=speed_acc, mode='lines', name='Vehicle Speed (km/h)', line=dict(color='blue', width=3)),
    secondary_y=False
)
fig3.add_trace(
    go.Scatter(x=time_acc, y=disp_acc, mode='lines', name='Displacement (m)', line=dict(color='red', width=2, dash='dash')),
    secondary_y=True
)

# Actual 50 km/h marker
idx_50 = np.argmax(speed_acc >= 50)
if idx_50 > 0:
    t_50 = time_acc[idx_50]
    fig3.add_vline(x=t_50, line_width=1, line_dash="dot", line_color="orange", opacity=0.7)
    fig3.add_annotation(x=t_50, y=50, text=f"Actual 50km/h @ {t_50:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=-30)

# Actual top speed marker
idx_max = np.argmax(speed_acc >= speed_kmh * 0.99)
if idx_max > 0:
    t_max = time_acc[idx_max]
    fig3.add_vline(x=t_max, line_width=1, line_dash="dot", line_color="green", opacity=0.7)
    fig3.add_annotation(x=t_max, y=speed_kmh, text=f"Actual Top Speed @ {t_max:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=30)

# Target 0-50 km/h line
fig3.add_vline(x=accel_time_0to50, line_width=1, line_dash="dot", line_color="purple", opacity=0.7)
fig3.add_annotation(x=accel_time_0to50, y=50, text=f"Target 50km/h @ {accel_time_0to50:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=-50, font=dict(color="purple"))

# Target 0-top speed line
fig3.add_vline(x=accel_time_full, line_width=1, line_dash="dot", line_color="brown", opacity=0.7)
fig3.add_annotation(x=accel_time_full, y=speed_kmh, text=f"Target Top Speed @ {accel_time_full:.1f}s", showarrow=True, arrowhead=2, ax=20, ay=-70, font=dict(color="brown"))

fig3.update_layout(
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
    margin=dict(l=20, r=20, t=40, b=20),
    height=400
)
fig3.update_xaxes(title_text="Time (s)")
fig3.update_yaxes(title_text="Speed (km/h)", secondary_y=False)
fig3.update_yaxes(title_text="Displacement (m)", secondary_y=True)
st.plotly_chart(fig3, use_container_width=True)

st.markdown("---")
st.caption("💡 Tip: Purple dashed line = target 0→50 km/h time, brown dashed line = target 0→top speed time.")