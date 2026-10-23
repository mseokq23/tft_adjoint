# 실험 로그 — 오픈소스 TCAD(DEVSIM) 산화물 TFT 보정 파라미터의 식별가능성·최적 측정조건

- 기준 논문: G. Jeong et al., arXiv:2609.12184 (agentic TCAD calibration, BG-IWO TFT, Sentaurus)
- 모델: 같은 descriptor 클래스(S/D·게이트 WF, Qf, ND,eff, Gaussian donor + 지수 tail acceptor DOS, µ0), 전자 DD, 노드 5,889개
- 기본 잡음 가정: 디지타이즈한 log10 ID의 σ = 0.02 dec(10 pA 이상만 사용). Fisher 결과는 이 가정과 Step-5 근방 선형화에 의존하고, σ에 정비례함. 실제 모델 오차는 이보다 크므로(논문 hold-out 오차 Vth 25–46 mV) 여기 결론은 낙관적(보수적) 쪽
- 실행 명령과 결과 폴더는 README.md 참고
- 최종 업데이트: 2026-10-23

## 합성 데이터(가상 실측) 결과

| 실험 | 핵심 결과 | 파일 |
|---|---|---|
| E0 기울기 검증 | adjoint vs 재계산 FD normwise 상대오차 ≤ 4.4e-6(아문턱), ≤ 3e-8(on/선형). Taylor 1.00/2.00. 바이어스점당(11개 θ) adjoint 0.72–0.76 s vs FD 2.9–5.8 s → 4–8배. `requirements.txt`로 새로 만든 환경에서 다시 돌려 시간 외 모든 값이 이전 실행과 같음을 확인 | results/e0 |
| E1 순방향 재현 | 논문 Step0→5, 4%-W, L = 2 µm 변화 방향 8/8 일치. Step5(VD = 1 V): Vth −0.177 V, SS 122.7 mV/dec, Ion 5.55e-4 A. 4%-W: Vth +0.498 V, Ion ×0.63. DOS 격자 ΔVth ≤ 0.7 mV, 메쉬 2배 ΔVth 0.07 mV·Ion −0.6 % | results/e1, e1_w4 |
| E3/E3b Fisher 식별가능성 | Qf ≡ 게이트 WF 정확 축퇴(고유값 비 2.2e-17; Qf 1e11 cm⁻² = WF 6.3 mV = q/C_ox). Qf 제외 시 조건수 8.0e7, 1σ: ND,eff ±2.95 dec, N_TA ±2.37 dec, N_GD ±1.92 dec, WF(eff) ±0.98 eV. DOS 고정 시 WF ±8.8 mV·ND ±0.044 dec | results/e3 |
| E3c 강건성 | 4%-W 세트, Nc = 2e18, ε_ch = 15에서도 Qf–WF 정확 축퇴와 ND,eff 중심 sloppy 방향 유지. 4%-W 세트에서는 ND,eff 민감도가 2%-W의 1/50(RMS 0.014 vs 0.66 dec/decade), 1σ 7,300 dec → 사실상 무관 | results/e3c |
| 논문 hold-out의 정보량 | VD 0.1/0.2/0.5 V·L = 2 µm 5개 곡선을 다 더해도 DOS/ND/WF σ는 0.47–0.60배(반복 측정 수준 0.53배). L = 2 µm는 wf_sd·µ0 분리에만 효과(0.26–0.30배) | results/e3/subspace.md |
| **E5 profile likelihood (비선형)** | ND,eff를 10^15.9–10^18.9 cm⁻³(1000배)로 고정하고 나머지 8개를 재보정해도 χ² ≤ 0.28(95 % 기준 3.84), RMS ≤ 0.0013 dec. 이 세트들의 논문 hold-out 5곡선 예측도 진짜와 RMS ≤ 0.0034 dec, ΔVth ≤ 0.8 mV(논문 자체 hold-out 오차 25–46 mV의 1/30 이하). 10^19.4(χ² 188)에서도 hold-out ΔVth는 ≤ 13.7 mV로 논문 검증 오차보다 작음. 10^19.15에서 χ² 4.2 → 95 % 상한 약 1.4e19, 하한 없음(10^15.9에서도 χ² 0.003; ND,eff가 작아지면 영향 자체가 사라짐). 2%-W 데이터는 4%-W 논문값(6.4e16)과 같은 ND,eff로도 똑같이 맞음 | results/e5_nd |
| **E8 트랩 종류 축퇴(해석적+수치)** | N(1−f) = N − N·f이므로 donor형 DOS ≡ 같은 모양의 acceptor형 DOS + 균일 고정전하(ND,eff에 흡수). Gaussian을 acceptor형으로 바꾸고 ND,eff를 4.1e18 올리면 I–V 7곡선 max Δlog ID 3.6e-15 dec, C–V 1.7e-14 → 기계 정밀도로 동일. ND,eff를 함께 보정하는 한 donor/acceptor 구분은 원리적으로 불가 | results/e8 |
| E4 공간·에너지 지도 | adjoint 1회(9–14 ms)로 채널 노드 3,344개 국소 민감도, 합 검증 1.000000. DOS 민감도는 페르미 준위 능선을 따라감 | results/e4map, e4energy |
| E6 DOS 검출 창 | 10 pA 하한에서 E_F(전면 계면) = Ec − 0.35 eV. DOS가 검출 한계보다 큰 구간 Ec − 0.05–0.35 eV(acceptor), 0.05–0.30 eV(donor). ND,eff를 함께 보정하면 두 종류의 검출 한계가 일치(E8과 같은 이유) | results/e6 |
| E4/E7 확장 OED | 2곡선 추가 기준 Δlog det F: t_IWO 7 nm 27.3, t_IWO 5 nm 22.7, t_ox 14 nm 22.4, 준정적 C–V(0.5 % 잡음) 21.0(2 %면 12.1), t_ox 10 nm 18.5 vs 논문 L = 2 µm 7.4, 논문 hold-out 5곡선 전체 12.4. Qf–WF 축퇴는 t_ox 변화로만 해소(σ_Qf 2.5e11, σ_WF 5.4 mV). t_ox14 + t_IWO7 + C–V: σ N_TA 2.37 → 0.14 dec, ND 2.95 → 0.44 dec, Qf 1.1e10 cm⁻². C–V 범함수는 DEVSIM 접촉전하와 1e-16, 기울기는 FD와 ≤ 2e-6 일치 | results/e4, e7 |
| **E9 측정 추천의 비선형 검증** | E5에서 보정·검증 데이터로 구분되지 않던 세트(ND,eff 10^15.9–10^18.9)를 새 측정으로 예측: t_IWO 7 nm I–V χ² 2.4e3–3.0e4(ΔVth 34–190 mV), t_ox 14 nm I–V χ² 523–7.5e3(ΔVth 18–71 mV), 준정적 C–V(0.5 % 잡음) χ² 51–1,380 → 모두 기각(기준 3.84). C–V는 2 % 잡음이면 가장 가까운 세트(10^17.9)를 못 가름(χ² 3.2) | results/e9 |
| 그림 | results/fig_abstract.png: (a) ND,eff profile + hold-out 예측, (b) 1000배 다른 ND,eff의 같은 곡선, (c) 측정 추가별 1σ, (d) DOS 검출 창 | scripts/fig_abstract.py |

## 실측(디지타이즈) 비교

| 단계 | 결과 | 파일 |
|---|---|---|
| 디지타이즈 | 논문 그림은 비트맵. 반투명 마커 색 + 측정 VG 격자(0.05 V)를 이용해 마커 중심을 읽는 디지타이저 작성(scripts/digitize_raster.py, digitize_paper.py). 9개 곡선 추출. 같은 곡선이 여러 그림에 나올 때 그림 간 차이 RMS 0.028 dec → 판독 σ ≈ 0.02 dec(분석에 가정한 값과 일치). VD 0.1/0.2/0.5와 4%-W VD 0.05의 아문턱은 마커가 겹쳐 점이 적음(신뢰도 낮음) | data/meas.csv, data/digitize/consistency.md |
| E10 논문 파라미터 → DEVSIM | Step 5 값 그대로: VD = 1 V ΔVth +208 mV, E_log 1.2 dec(논문 Sentaurus 31.5 mV). 단계별 J도 단조 감소하지 않음(Step 4에서 악화). descriptor 값은 시뮬레이터·재료 상수(χ 등)에 묶인 유효값 | results/e10 |
| E11 실측 보정(2%-W) | 우리 모델을 실측 VD 0.05/1 V에 adjoint-LM 보정(Qf 고정, 바닥 잡음 구간 제외): RMS 0.11 dec(σ_eff 0.12 dec, 모델 오차가 판독 오차의 6배), ΔVth −2/+40 mV, ΔSS +23/+21 mV/dec, ΔIon ±0.017 dec. 남는 차이는 아문턱 VD 의존성(실측 ~100 mV, 모델 ~60 mV). 실측 FIM 1σ: ND,eff·DOS 수백 dec(정보 사실상 없음) | results/e11_2W_v2 |
| E11 hold-out(실측) | 보정 모델 예측: VD 0.1/0.2/0.5 ΔVth −34/−38/+58 mV, L = 2 µm VD 0.05/1 +4/+90 mV, ΔIon ≤ 0.016(L1)·0.08–0.11 dec(L2). 논문(Sentaurus) 24.9/36.1 mV, 27.4/46.2 mV보다 다소 큼 | results/e11_2W_v2 |
| **E12 실측 profile** | 실측 2%-W 보정 곡선으로 ND,eff를 10^16.85–10^19.05(약 160배)에 고정하고 나머지 8개 재보정: Δχ²_eff ≤ 0.77(95 % 기준 3.84, σ_eff 0.12 dec 기준). 논문 4%-W 값(6.4e16)과 같은 수준의 ND,eff로도 2%-W 실측이 똑같이 맞음. 이 세트들의 실측 hold-out ΔVth는 조건별로 거의 같음(VD 0.1: −31…−34, VD 0.2: −38…−39, L2 VD 0.05: +3…+9, L2 VD 1: +90…+94 mV; 점이 적은 VD 0.5만 +48…+70 mV) → 실측 hold-out으로도 구분 불가. 10^19.2 이상은 수치 수렴 실패로 미계산 | results/e12_nd, results/fig_measured.png |
| 4%-W | 4%-W VD 0.05는 문턱 구간(−0.25…0.45 V) 판독 불가, VD 1은 바닥 잡음 혹 → 우리 모델 보정 RMS 0.33 dec로 불량. 4%-W 실측 결론은 보류(대신 E12가 2%-W 데이터만으로 4%-W 수준 ND,eff를 허용함을 보임) | results/e11_4W_v2 |

## 해석할 때 주의할 점

- 결론은 "같은 descriptor 클래스"와 이 오픈소스 구현(DEVSIM, 전자 DD, 가정한 재료 상수)에 대한 것입니다. 기준 논문의 시뮬레이터(Sentaurus)와는 같은 파라미터에서 Vth가 약 0.2 V 다릅니다(E10).
- 실측 비교에서 모델 오차(σ_eff 0.12 dec)가 판독 오차(0.02 dec)보다 크고, 아문턱 VD 의존성을 과소 예측합니다(E11).
- 4%-W 곡선과 VD 0.1/0.2/0.5 V 아문턱은 판독 한계로 보조 자료입니다.
- adjoint TCAD 자체는 새 방법이 아닙니다(DEVSIM 2009, Charon/Sacado, Burger–Pinnau 2003, Khalil 1995).
