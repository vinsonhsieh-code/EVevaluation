# ================== 圖2：車輪推力 vs 車速 + 工作點 ==================
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

# ---------- 新增：平路負載線與最大車輪推力曲線的交點 ----------
# 注意：v_from_n 對應的 x 是車速 (km/h)，F_wheel_max 對應 y；speed_kmh_flat 對應 x，F_wheel_flat 對應 y
intersections_flat_force = find_intersection(v_from_n, F_wheel_max, speed_kmh_flat, F_wheel_flat)
for x_cross, y_cross in intersections_flat_force:
    fig2.add_trace(go.Scatter(
        x=[x_cross], y=[y_cross], mode='markers',
        marker=dict(color='red', size=10, symbol='x'),
        name='平路交點', showlegend=False,
        text=f"平路交點<br>車速: {x_cross:.1f} km/h<br>推力: {y_cross:.0f} N",
        hoverinfo='text'
    ))

# ---------- 新增：爬坡負載線與最大車輪推力曲線的交點（若有） ----------
if F_wheel_climb is not None:
    intersections_climb_force = find_intersection(v_from_n, F_wheel_max, speed_kmh_climb, F_wheel_climb)
    for x_cross, y_cross in intersections_climb_force:
        fig2.add_trace(go.Scatter(
            x=[x_cross], y=[y_cross], mode='markers',
            marker=dict(color='green', size=10, symbol='x'),
            name='爬坡交點', showlegend=False,
            text=f"爬坡交點<br>車速: {x_cross:.1f} km/h<br>推力: {y_cross:.0f} N",
            hoverinfo='text'
        ))

x_max = max(v_max_motor, speed_kmh) * 1.15 if max(v_max_motor, speed_kmh) > 0 else 100
fig2.update_yaxes(title_text="車輪推力 (N)", range=[y_min_force, y_max_force], tickvals=y_ticks_f, tickfont=dict(color='white'), zeroline=True, zerolinecolor='gray', zerolinewidth=1.5)
fig2.update_xaxes(title_text="車速 (km/h)", range=[0, x_max], tickfont=dict(color='white'), zeroline=True, zerolinecolor='gray', zerolinewidth=1.5)

fig2.add_annotation(x=0, y=F_wheel_peak, xref="x", yref="y", text=f"<b>{F_wheel_peak:.0f}</b>", showarrow=False, xanchor="right", xshift=-15, font=dict(color="dodgerblue", size=14), bgcolor="rgba(26,28,35,0.9)", bordercolor="dodgerblue", borderwidth=1, borderpad=4)
fig2.add_annotation(x=speed_kmh, y=F_design_flat, xref="x", yref="y", text=f"<b>目標: {speed_kmh:.0f} km/h</b>", showarrow=True, arrowhead=2, arrowcolor="orange", arrowsize=1, arrowwidth=2, ax=-55, ay=-65, font=dict(color="orange", size=12), bgcolor="rgba(26,28,35,0.9)", bordercolor="orange", borderwidth=1, borderpad=3)
fig2.add_annotation(x=v_max_motor, y=F_at_vmax, xref="x", yref="y", text=f"<b>極速: {v_max_motor:.0f} km/h<br>{F_at_vmax:.0f} N</b>", showarrow=True, arrowhead=2, arrowcolor="purple", arrowsize=1, arrowwidth=2, ax=65, ay=-45, font=dict(color="#d8b4e2", size=12), bgcolor="rgba(26,28,35,0.9)", bordercolor="purple", borderwidth=1, borderpad=3)

fig2.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5), margin=dict(l=80, r=110, t=100, b=20), height=550)
st.plotly_chart(fig2, use_container_width=True)
st.markdown("---")