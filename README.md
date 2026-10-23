# tft_adjoint — 오픈소스 TCAD(DEVSIM) 산화물 TFT 보정·식별가능성 연구 코드

KCS 2027 초록 주제: **"Adjoint 민감도 기반 오픈소스 TCAD를 이용한 산화물 반도체 TFT 보정, 파라미터 식별가능성 분석 및 최적 측정조건 설계"**

비교 기준 논문: G. Jeong et al., *Agentic TCAD Calibration Workflow for Oxide Semiconductor Transistors*, arXiv:2609.12184 (Georgia Tech·Samsung·NVIDIA). 이 논문의 BG-IWO TFT 구조와 보정 파라미터(일함수, Qf, ND,eff, Gaussian·tail DOS, µ0)를 같은 모델 클래스로 DEVSIM에 구현했습니다.

## 설치

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt        # 검증 환경: Python 3.13.16, Linux x86-64
```

- DEVSIM은 import할 때 BLAS/LAPACK 라이브러리를 찾습니다. `requirements.txt`의 `mkl`(Intel MKL, x86-64 전용)이 그 역할을 합니다. MKL을 쓸 수 없는 환경에서는 BLAS/LAPACK 라이브러리를 설치하고 환경변수 `DEVSIM_MATH_LIBS`로 지정해야 합니다(DEVSIM 문서 참고). 검증은 MKL로만 했습니다.
- DEVSIM 기본 직접 솔버(MKL PARDISO)는 일부 예제에서 실패해, DEVSIM 휠에 들어 있는 UMFPACK shim(`devsim.umfpack`)을 기본으로 씁니다(`TFT(..., solver="umfpack")`). 따로 설치할 것은 없습니다.
- 스크립트는 모두 레포 루트에서 실행합니다(입력 파일 경로가 루트 기준).

## 데이터 출처

- `data/meas.csv`: arXiv:2609.12184 (v2) Fig. 5(b), 6(a), 6(b), 7(a), 7(c)에 실린 측정 전달곡선을 `scripts/digitize_paper.py`로 디지타이즈한 값입니다.
  열은 `set`(2W/4W), `L_um`, `vd`(V), `vg`(V), `id`(A), `n_panels`(그 점을 읽은 그림 수), `spread_dec`(그림 간 차이, dec)입니다.
  같은 곡선이 여러 그림에 나오면 중앙값을 썼고, 그림 간 일치도는 `data/digitize/consistency.md`에 있습니다(판독 σ ≈ 0.02 dec).
- `data/meas_panels.csv`: 그림(패널)별 판독값.
- 측정 데이터의 출처는 위 논문입니다. 사용할 때는 원 논문을 인용해 주세요.
- 논문 PDF, 그 안의 그림 비트맵, 비트맵 위에 판독점을 겹친 확인용 그림은 저작권 때문에 넣지 않았습니다.
  디지타이즈를 다시 하거나 확인하려면 arXiv에서 v2 PDF를 받아 `data/refs/jeong2026_agentic_tcad.pdf`로 저장한 뒤 아래를 실행합니다.
  같은 판의 PDF라면 `data/meas.csv`, `data/meas_panels.csv`, `data/digitize/consistency.md`가 바이트 단위로 똑같이 다시 만들어집니다.

```bash
python scripts/extract_figs.py       # → data/figs/p5_img0.png, p5_img1.png, p6_img0.png (MD5로 같은 판인지 확인)
python scripts/digitize_paper.py     # → data/meas.csv, data/meas_panels.csv, data/digitize/ (약 30초)
python scripts/digitize_check.py     # 판독점을 원 그림 위에 겹쳐 그림 → data/digitize/check_*.png
```

## 구성

**`tftsim/` — 모델·민감도 라이브러리**

| 파일 | 내용 |
|---|---|
| `tftsim/params.py` | 구조(Table III), 재료 상수(가정값), 논문 Step 0–5·4%-W 파라미터, θ 정의, DOS 이산화 |
| `tftsim/device.py` | 2D 메쉬, 전자 drift-diffusion, 이산화 트랩, ND,eff, Qf, 일함수 접촉, 바이어스 램프·스윕 |
| `tftsim/sensitivity.py` | Jacobian 추출, 잔차 차분 ∂F/∂θ, adjoint·forward 민감도(스케일링 LU + 반복 보정), 선형 범함수 기울기, 재계산 FD 기준값 |
| `tftsim/metrics.py` | Vth(100 nA), SS(0.1–100 nA), Ion 추출과 그 기울기, 목적함수 J(자체 정규화 계수) |
| `tftsim/calib.py` | 바이어스점별 warm start로 곡선을 다시 계산하는 CurveModel, 로그 전류 잔차 문제 |
| `tftsim/cv.py` | 게이트 전하 Q_G(전위의 선형 범함수)와 adjoint 기울기, 준정적 C = dQ_G/dV_G |

**합성 데이터(가상 실측) 실험 — 진짜 = 논문 Step-5 파라미터**

| 스크립트 | 내용 | 결과 |
|---|---|---|
| `scripts/e0_gradient_check.py` | E0 기울기 검증: adjoint·forward·재계산 FD 비교, Taylor 테스트, FD 스텝 스캔 | `results/e0` |
| `scripts/e1_forward_reproduce.py` | E1 순방향 재현: Step 0–5, 4%-W, L = 2 µm, 수치 수렴 점검 | `results/e1`, `e1_w4` |
| `scripts/e2_twin_calibration.py` | E2 가상 실측 보정: adjoint-LM 대 유한차분 LM(대 DFO-LS) | 결과 미포함 |
| `scripts/e3_identifiability_oed_preview.py` | E3 Fisher 정보 기반 식별가능성, 논문 검증(hold-out) 조건의 OED 순위 | `results/e3` |
| `scripts/e3b_identifiable_subspace.py` | E3b 정확히 축퇴된 Qf 방향을 뺀 부분공간에서 식별가능성·OED 재계산(시뮬레이션 없음) | `results/e3/subspace.*` |
| `scripts/e3c_robustness.py` | E3c 4%-W 세트·다른 재료 가정(Nc, ε)에서도 식별가능성 구조가 유지되는지 | `results/e3c` |
| `scripts/e4_sensitivity_map.py` | E4 adjoint 1회로 얻는 공간 민감도 지도(채널 어디의 전하가 ID를 바꾸는가) | `results/e4map` |
| `scripts/e4_energy_sensitivity.py` | E4 에너지별 트랩 민감도 지도(전달곡선이 DOS의 어느 에너지를 보는가) | `results/e4energy` |
| `scripts/e4_oed_extra_conditions.py` | E4 논문 밖 측정 조건(게이트 산화막 10·14 nm)이 Qf–게이트 일함수 축퇴를 푸는지 | `results/e4` |
| `scripts/e5_profile_likelihood.py` | E5 비선형 profile likelihood(ND,eff 고정 후 나머지 재보정) + 재보정 세트의 hold-out 예측. `--init`으로 이전 결과에서 이어 달리기 | `results/e5_nd_*` |
| `scripts/e5b_merge_profiles.py` | 나눠 돌린 profile 결과 병합(표·그림) | `results/e5_nd` |
| `scripts/e6_dos_detectability.py` | E6 에너지 구간별 DOS 검출 한계(교란 파라미터 profiling 포함)와 페르미 준위 창 | `results/e6` |
| `scripts/e7_oed_cv_tch.py` | E7 준정적 C–V, 채널 두께 변형까지 포함한 확장 OED | `results/e7` |
| `scripts/e8_trap_type_equivalence.py` | E8 donor형 DOS ≡ acceptor형 DOS + 균일 고정전하(ND,eff) 수치 증명(I–V, C–V) | `results/e8` |
| `scripts/e9_discrimination_test.py` | E9 E5에서 구분되지 않던 세트들을 추천 측정(t_IWO 7 nm, t_ox 14 nm, C–V)이 기각하는지 비선형 확인 | `results/e9`, `e9a`, `e9b` |
| `scripts/fig_abstract.py` | 초록용 4패널 그림(결과 파일만 읽음) | `results/fig_abstract.png` |

**실측(디지타이즈) 비교**

| 스크립트 | 내용 | 결과 |
|---|---|---|
| `scripts/extract_figs.py` | 논문 PDF에서 그림 비트맵 추출(PDF는 직접 받아야 함) | `data/figs/`(미포함) |
| `scripts/digitize_raster.py`, `scripts/digitize_paper.py` | 비트맵에서 실측 곡선 추출(반투명 마커 색 + 0.05 V 측정 격자), 그림 간 일치도 | `data/meas.csv`, `data/meas_panels.csv`, `data/digitize/consistency.md` |
| `scripts/digitize_check.py` | 추출한 곡선을 원래 그림 위에 다시 그려 육안 확인 | `data/digitize/check_*.png`(미포함) |
| `scripts/e10_measured_comparison.py` | E10 논문 파라미터를 DEVSIM에 넣었을 때 실측과의 차이(J, Vth, SS, Ion) | `results/e10` |
| `scripts/e11_calibrate_measured.py` | E11 실측 곡선에 우리 모델을 adjoint-LM으로 보정, hold-out 실측 예측 | `results/e11_*` |
| `scripts/e12_profile_measured.py` | E12 실측 데이터 기준 ND,eff profile likelihood | `results/e12_nd` |
| `scripts/fig_measured.py` | 실측 비교 4패널 그림 | `results/fig_measured.png` |
| `scripts/digitize_utils.py` | 벡터 그림 경로·WebPlotDigitizer CSV를 meas.csv 형식으로 바꾸는 초기 도구. 이 논문 그림은 비트맵이라 최종 데이터에는 쓰지 않음 | – |

`scripts/_common.py`는 공용 도우미(경로, JSON 저장, 타이머)이고, `EXPERIMENT_LOG.md`에 실험별 핵심 수치가 있습니다.

## 실행 (results/의 결과 재현)

시간은 2 vCPU 기준 실측값입니다. 결과는 `results/<실험>/`에 `summary.md`, `summary.json`, CSV/NPZ, PNG로 저장됩니다.

```bash
# 합성 데이터 실험
python scripts/e0_gradient_check.py                 # 약 2분
python scripts/e1_forward_reproduce.py              # 약 30분 (--quick: 약 5분)
python scripts/e1_forward_reproduce.py --only w4 --out results/e1_w4   # 약 4분
python scripts/e3_identifiability_oed_preview.py    # 약 10분
python scripts/e3b_identifiable_subspace.py         # 수 초 (E3 결과 후처리)
python scripts/e3c_robustness.py                    # 약 6분
python scripts/e4_sensitivity_map.py                # 약 1분
python scripts/e4_energy_sensitivity.py             # 약 2분
python scripts/e4_oed_extra_conditions.py           # 약 9분
python scripts/e6_dos_detectability.py              # 약 2분
python scripts/e7_oed_cv_tch.py                     # 약 10분 (E3·E4 결과 필요; --reuse: 저장된 행 재사용)
python scripts/e8_trap_type_equivalence.py          # 약 10분
# E5: 진짜 값(10^18.41)에서 양쪽으로 continuation. 두 쪽을 따로 돌린 뒤 병합 (각 1–2시간)
python scripts/e5_profile_likelihood.py --values 18.9,19.15,19.4 --out results/e5_nd_up
python scripts/e5_profile_likelihood.py --values 17.9,17.4,16.9,15.9 --out results/e5_nd_down
python scripts/e5b_merge_profiles.py results/e5_nd_up results/e5_nd_down --out results/e5_nd
# E9: E5 세트 중 χ² < 3.84인 것만 사용. 둘로 나눠 돌린 뒤 병합
python scripts/e9_discrimination_test.py results/e5_nd --values 15.9,16.9,17.4 --out results/e9a
python scripts/e9_discrimination_test.py results/e5_nd --values 17.9,18.9 --out results/e9b
python scripts/e9_discrimination_test.py results/e5_nd --combine results/e9a results/e9b --out results/e9
python scripts/fig_abstract.py

# 실측 비교 (data/meas.csv 포함. 디지타이즈부터 다시 하려면 '데이터 출처' 참고)
python scripts/e10_measured_comparison.py           # 약 7분
python scripts/e11_calibrate_measured.py --set 2W --start step5 --out results/e11_2W                          # 1단계, 약 15분
python scripts/e11_calibrate_measured.py --set 2W --start results/e11_2W/summary.json --out results/e11_2W_v2   # 2단계(보고 값), 약 30분
python scripts/e11_calibrate_measured.py --set 4W --vgmin -0.7 --start w4 --out results/e11_4W                 # 약 35분
python scripts/e11_calibrate_measured.py --set 4W --vgmin -0.7 --start results/start_4W.json --out results/e11_4W_v2
python scripts/e12_profile_measured.py --best results/e11_2W_v2/summary.json --values 18.35,17.85,17.35,16.85 --out results/e12_nd   # 점당 약 15–20분
python scripts/e12_profile_measured.py --best results/e11_2W_v2/summary.json --values 19.05 --out results/e12_nd
python scripts/fig_measured.py
```

참고:

- 저장된 E5 결과는 같은 점들을 네 번에 나눠(`e5_nd_up`, `_up2`, `_down`, `_down2`; 일부는 `--init`으로 이전 실행에서 이어서) 계산한 뒤 `e5b_merge_profiles.py`로 네 폴더를 합친 것입니다. 위처럼 한 번에 돌리면 continuation 경로가 달라 마지막 자릿수는 다를 수 있습니다.
- E11 2%-W 1단계(`results/e11_2W`)는 바닥 잡음 구간 처리(`--floor`) 규칙을 넣기 전에 돌렸습니다. 지금 코드로 다시 돌리면 1단계 수치는 달라질 수 있습니다. 보고 값은 2단계(`results/e11_2W_v2`)입니다.
- `results/start_4W.json`은 논문 4%-W 값에서 게이트 일함수만 4.45 → 4.25 eV로 바꾼 시작점입니다(E10에서 본 약 +0.2 V Vth 차이만큼).
- E12에서 ND,eff = 10^19.2는 수치 수렴에 실패해 10^19.05로 대신 계산했습니다. `e12_profile_measured.py`는 실행할 때마다 `results/e12_nd/` 아래의 모든 점을 모아 표와 그림을 다시 만듭니다.
- E2는 코드만 있고 결과는 넣지 않았습니다.
- 실행 로그(`results/*_log.txt`)는 넣지 않았습니다.
- 동시에 두 스크립트를 돌려도 되지만, 한 프로세스 안에서 TFT를 새로 build하면 DEVSIM 전역 상태가 초기화되므로(`reset_devsim`) 한 프로세스에서는 장치를 하나씩만 다룹니다.

## 모델 가정과 한계 (논문에 없는 값은 가정)

- **재료 상수(논문 미공개):** χ = 4.16 eV, Eg = 3.05 eV, Nc = 5×10¹⁸ cm⁻³, ε(IWO) = 10 (Fung et al., JAP 2009 a-IGZO 값), ε(HfO2) = 20.
  모델에서는 일함수가 항상 (WF − χ) 형태로만 들어가므로 **χ와 일함수는 완전히 축퇴**됩니다. 논문의 절대 일함수 값은 우리 χ 가정 아래에서만 의미가 있습니다.
- 전자만 고려하는 drift-diffusion(Scharfetter–Gummel), Boltzmann 통계, 상수 이동도.
- 트랩 점유는 전도대와 준평형(f = n/(n+n1))이라고 가정했습니다. DC 스윕에서만 유효하고 과도·히스테리시스는 다루지 않습니다.
- DOS는 Ec 아래 0–1.6 eV를 40개 준위로 이산화했습니다(0.025 eV 간격 32개 + 0.1 eV 간격 8개). Gaussian은 Ec에서 잘랐습니다(`trap_levels`).
- S/D: 평형 경계조건 Schottky 접촉(n = Nc·exp(−(WF−χ)/Vt), 무한 표면 재결합 속도). 게이트: HfO2 위 Dirichlet.
- 패시베이션 HfO2는 채널 위 간격 영역만 모델링했고, S/D 금속 측벽과 패시베이션의 접촉은 무시했습니다.
- 기울기 ∂F/∂θ는 수렴해에서 잔차 조립 1–2회를 차분해 얻습니다(재계산 없음). 기호 미분이 아니므로 "exact"라고 부르지 않습니다.

## 수치 설정에서 알아둘 것

- DEVSIM `solve`는 절대오차와 상대오차를 **모두** 만족해야 수렴합니다. 여기서는 absolute_error = 1e30으로 두어 상대오차 1e-10만 사용합니다.
- 아문턱 영역에서는 상대오차가 1e-9 근처에서 몇 회 정체한 뒤 수렴합니다. 그래서 최대 반복을 100회로 두었습니다. 확장 정밀도(float128)로도 정체는 바뀌지 않았습니다. 초기 build가 이 기준에서 수렴하지 않으면 완화한 기준(상대오차 1e-8, 최대 200회)으로 한 번 더 시도합니다.
- 민감도 선형계는 전위(V)와 전자 농도(최대 ~1e20 cm⁻³)가 섞여 조건이 나쁩니다. 일반 LU를 쓰면 출력 성분에서 5–6자리를 잃으므로, 행·열 equilibration과 반복 보정 2회(`ScaledLU`)를 적용했습니다.
