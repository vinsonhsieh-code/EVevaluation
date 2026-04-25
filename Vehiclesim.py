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

        # ---------- 追加實際計算範例（動態帶入目前參數與結果）----------
        st.markdown("---")
        st.markdown("### 📊 本次工況實際計算範例")
        
        # 取得目前側邊欄與計算結果中的關鍵數值
        try:
            # 基本參數
            m_kg = total_mass
            fr_val = fr
            cd_val = cd
            area_val = area
            gear_ratio_val = gear_ratio
            gear_eff_val = gear_eff / 100.0
            battery_kwh = user_battery_energy_kwh
            soc = battery_soc
            
            # 能耗計算結果（已在上方計算）
            dist_km = total_dist_km
            wh_km_batt = wh_per_km_batt
            drive_wh = drive_energy_wh
            regen_wh = regen_energy_wh
            net_wh = total_energy_wh
            range_km = estimated_range_km
            
            # 馬達峰值扭矩（若有使用效率地圖）
            if 'peak_torque_input' in locals() and peak_torque_input is not None:
                motor_torque_peak = peak_torque_input
            elif 'manual_peak_torque' in locals() and manual_peak_torque is not None:
                motor_torque_peak = manual_peak_torque
            else:
                motor_torque_peak = T_peak  # 來自馬達規格
            
            # 效率地圖使用狀態
            use_eff_map = (eff_interpolator is not None)
            
            st.markdown(f"""
            **🔧 參數盤點與假設**  
            - **總質量**: {m_kg} kg (車重 {weight} kg + 載重 {load} kg)  
            - **阻力設定**: $C_d = {cd_val}$, 迎風面積 $A = {area_val} \, \\text{{m}}^2$, 滾阻 $f_r = {fr_val}$  
            - **傳動幾何**: 輪胎半徑 {wheel_radius_m:.4f} m, 齒輪比 {gear_ratio_val:.1f}, 齒輪效率 {gear_eff}%  
            """)
            
            if use_eff_map:
                st.markdown(f"""
            - **馬達效率地圖解析**: 將上傳的效率地圖 CSV 中「{df_eff_converted.columns[0]}」轉速欄位與各扭矩百分比欄位（10%~100%）自動對應到您提供的**峰值扭矩 {motor_torque_peak} Nm** 的百分比（即 Y 軸: {motor_torque_peak*0.1:.1f} Nm ~ {motor_torque_peak} Nm），建立 2D 效率插值網格。  
                """)
            else:
                st.markdown(f"""
            - **馬達效率**: 使用固定效率 {motor_eff}% (未上傳效率地圖)。  
                """)
            
            st.markdown(f"""
            - **可用電池能量**: {battery_kwh} kWh × {soc}% = {battery_kwh * soc / 100:.2f} kWh (即 {battery_kwh * soc / 100 * 1000:.0f} Wh)  

            **📈 [IDC 工況逐秒積分運算結果]**  
            - **工況單圈總距離**: {dist_km:.3f} km  
            - **驅動耗電 (電池端)**: {drive_wh / dist_km:.2f} Wh/km (總計 {drive_wh:.2f} Wh)  
            - **回收充電 (電池端)**: {regen_wh / dist_km:.2f} Wh/km (總計 {regen_wh:.2f} Wh)  
            - **電池端單位淨能耗**: {wh_km_batt:.2f} Wh/km  
            - **理論預估續航里程**: {range_km:.1f} 公里  

            **✍️ 計算過程舉例（選取某個時間點，例如最高功率點）**  
            """)
            
            # 嘗試找出一個典型驅動點（功率較大的點）來說明效率查表與功率換算
            if 'times' in locals() and len(times) > 0 and 'power_batt_w' in locals():
                # 找驅動功率最大的時間點索引（正功率）
                pos_power_mask = power_batt_w > 0
                if np.any(pos_power_mask):
                    idx_max = np.argmax(power_batt_w[pos_power_mask])
                    orig_idx = np.where(pos_power_mask)[0][idx_max]
                    t_ex = times[orig_idx]
                    v_ex = df_energy['speed_kmh'].iloc[orig_idx]
                    a_ex = df_energy['accel_ms2'].iloc[orig_idx]
                    # 計算該點阻力
                    v_ms = v_ex / 3.6
                    F_roll = m_kg * G * fr_val
                    F_aero = 0.5 * RHO * cd_val * area_val * v_ms**2
                    F_acc = m_kg * a_ex
                    F_total = F_roll + F_aero + F_acc
                    P_wheel_ex = F_total * v_ms
                    # 取得馬達轉速與扭矩（來自工作點）
                    if df_operating is not None:
                        rpm_ex = df_operating['motor_rpm'].iloc[orig_idx]
                        torque_ex = df_operating['motor_torque_Nm'].iloc[orig_idx]
                    else:
                        rpm_ex = v_ms * 60 / (2 * math.pi * wheel_radius_m) * gear_ratio_val
                        torque_ex = F_total * wheel_radius_m / (gear_ratio_val * gear_eff_val)
                    
                    # 效率值
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
                    - 總效率：$\\eta_\\text{{total}} = {gear_eff_val*100:.1f}\\% \\times {eff_motor:.1f}\\% = {gear_eff_val * (eff_motor/100) * 100:.1f}\\%$  
                    - 電池端功率：$P_\\text{{batt}} = {P_wheel_ex:.0f} / {gear_eff_val * (eff_motor/100):.3f} = {power_batt_w[orig_idx]:.0f}\\,\\text{{W}}$  
                    - 該點瞬間能耗率即為 $P_\\text{{batt}}$，對全工況時間積分後得到總耗電 {net_wh:.1f} Wh。
                    """)
                    else:
                        st.markdown(f"""
                    在時間 **t = {t_ex:.1f} s** 時，車速 {v_ex:.1f} km/h，加速度 {a_ex:.2f} m/s²。  
                    - 行駛阻力：$F_\\text{{roll}} = {F_roll:.1f}\\,\\text{{N}}$，$F_\\text{{aero}} = {F_aero:.1f}\\,\\text{{N}}$，$F_\\text{{acc}} = {F_acc:.1f}\\,\\text{{N}}$  
                    - 總阻力：$F_\\text{{total}} = {F_total:.1f}\\,\\text{{N}}$  
                    - 輪上功率：$P_\\text{{wheel}} = {P_wheel_ex:.0f}\\,\\text{{W}}$  
                    - 固定馬達效率 {motor_eff}%，總效率 $\\eta_\\text{{total}} = {gear_eff}% \\times {motor_eff}% = {gear_eff_val * (motor_eff/100) * 100:.1f}\\%$  
                    - 電池端功率：$P_\\text{{batt}} = {P_wheel_ex:.0f} / {(gear_eff_val * (motor_eff/100)):.3f} = {power_batt_w[orig_idx]:.0f}\\,\\text{{W}}$  
                    - 該點瞬間能耗率即為 $P_\\text{{batt}}$，對全工況時間積分後得到總耗電 {net_wh:.1f} Wh。
                    """)
                    
                else:
                    st.markdown("> 無法找到典型的驅動功率點，跳過詳細計算舉例。")
            else:
                st.markdown("> 由於缺乏逐秒數據，無法展示詳細的時間點計算過程。")
                
        except Exception as e:
            st.markdown(f"> ⚠️ 動態計算範例時發生錯誤：{e}，請確認所有參數已正確設定。")