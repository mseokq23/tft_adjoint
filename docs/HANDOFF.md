# 인계 문서 — 새 세션에서 이어서 작업하기

이 문서는 작업을 새 환경이나 새 세션에서 정보 손실 없이 이어 가기 위한 것입니다. 연구 내용 해설은 [README.md](../README.md), 실험별 수치는 [EXPERIMENT_LOG.md](../EXPERIMENT_LOG.md)에 있습니다.

- 저장소: https://github.com/mseokq23/tft_adjoint (공개, 기본 브랜치 `main`)
- 기준 시점: 2026-10-25. 초록 마감: 2026-11-07(KCS 2027).
- 현재 상태: 1단계 가설검증 실험이 모두 끝났습니다(E0–E21: 결합 보정 E18, 일반 규칙 E19, 식별가능성 표 E20, 강건성 E21 포함). 가치 판단은 C4 점검(2026-10-24, 부분 겹침만 있음)까지 마쳐 초안 작성 조건을 충족했습니다(README 5장·6.3). 1쪽 초록 초안과 초록 그림, 논문용 자료(`docs/paper/`)를 만들었습니다. 2026-10-25에 기준 논문 저자의 공개 데모에서 같은 결론을 확인했고(E22), 초록 방향을 기준 논문 워크플로에 더하는 식별가능성 점검 단계로 정해 초록 그림을 바꿨습니다(v3). 남은 일은 초록 마무리(저자·소속 기입, 제출), 마감 전 선행연구 재확인, 전체 논문 작성, 다른 시뮬레이터 교차 확인입니다(README 6장).

## 1. 먼저 읽을 것 (순서대로)

1. README.md 1장(요약), 3장(결과 상세), 6장(앞으로 할 일)
2. `results/e14/summary.md`: 핵심 주장의 판정표(전달 시험, 퍼짐, 선형 기준, 강건성)
3. `results/e13/summary.md`, `results/e12_nd/summary.md`: 실측 등가 세트와 판별 예측
4. `results/crosscheck/README.txt`: 교차 확인용 세트 정의와 개념 설명
5. EXPERIMENT_LOG.md "현재까지의 평가"

## 2. 저장소 규칙

- **공개 저장소입니다.** 다음은 넣지 않습니다(`.gitignore`로 제외):
  - 논문 PDF(`data/refs/`)
  - 논문 그림 비트맵(`data/figs/`)
  - 비트맵에 판독점을 겹친 확인용 그림(`data/digitize/*.png`)
  - 실행 로그(`results/*_log.txt`)
- **커밋 작성자·커미터:** `mseokq23 <156693246+mseokq23@users.noreply.github.com>`. 저장소 로컬 설정으로 지정합니다.

  ```bash
  git config user.name mseokq23
  git config user.email 156693246+mseokq23@users.noreply.github.com
  git config commit.gpgsign false
  ```

- **커밋 메시지:** 내용만 씁니다. 공동 작성자(Co-Authored-By) 줄이나 자동으로 붙는 서명·링크 문구는 넣지 않습니다.
- **파일 내용:** 작업 환경 정보(편집·생성에 쓴 도구의 이름, 절대 경로 `/home/...`, `/tmp/...`)를 파일에 남기지 않습니다.
- **비공개 정보:** 연락 내용, 개인 연구 계획 메모는 저장소에 두지 않고 비공개 노트에 둡니다. 저장소에는 데이터 출처를 중립적으로만 씁니다.
- **push 전 점검:**

  ```bash
  git grep -n -i -E "co-authored-by|generated with|/home/|/tmp/" -- ':!docs/HANDOFF.md'   # 결과가 없어야 함
  git log --format='%an <%ae> | %cn <%ce>' | sort -u                 # 작성자·커미터가 위 한 가지뿐이어야 함
  ```

  편집·생성 도구의 이름이 파일, 커밋 메시지, PNG/NPZ 메타데이터에 들어가지 않았는지도 확인합니다.

## 3. 환경 복원

```bash
git clone https://github.com/mseokq23/tft_adjoint && cd tft_adjoint
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt     # Python 3.13, Linux x86-64에서 검증. mkl이 DEVSIM의 BLAS/LAPACK 역할
python scripts/e14_transfer_test_audit.py   # 수 초. 저장된 결과로 판정표 재생성(설치 확인용)
python scripts/e0_gradient_check.py         # 약 2분. DEVSIM과 민감도 엔진 확인
```

- `results/`에 지금까지의 결과가 모두 있습니다. 이어서 작업하려고 시뮬레이션을 다시 돌릴 필요는 없습니다.
- 시뮬레이션 없이 다시 만들 수 있는 결과 목록은 README 10장에 있습니다. 모두 저장된 파일과 바이트 단위로 같게 재생성됨을 확인했습니다(2026-10-24).
- 지금까지의 실행 환경은 2 vCPU였고, 백그라운드 작업 하나가 2시간을 넘지 않게 나눠 돌렸습니다.
  - profile 계산은 점 단위로 이어 달리기(`e5 --init`, e12의 continuation)가 되므로 나눠 실행할 수 있습니다.
  - 동시에 두 프로세스까지 돌렸습니다.

## 4. 핵심 수치와 출처

| 주장 | 수치 | 파일 |
|---|---|---|
| Qf ≡ 게이트 WF | Fisher 고유값 비 2.2e-17, Qf 1e11 cm⁻² = WF 6.3 mV | results/e3/summary.md, results/e4/summary.md |
| donor ≡ acceptor + ND,eff | I–V 최대 3.6e-15 dec, C–V 최대 1.7e-14(두께 변형 최대 8.1e-15) | results/e8, results/e8_thickness |
| 쌍둥이 비식별 범위 | ND,eff 10^15.9–10^18.9에서 χ² 최대 0.284, 95 % 상한 약 1.4e19 | results/e5_nd/summary.md |
| 쌍둥이 전달 시험 | max ΔVth ≤ 0.76 mV(10^19.4 세트도 13.7 mV) | results/e14/summary.md 1절 |
| 선형 기준 | 전달 시험 σ_pred/σ 0.25–0.43, 추천 시험 32–239 | results/e14/summary.md 3절 |
| 실측 비식별 범위 | ND,eff 10^16.85–10^19.05, Δχ²_eff 최대 1.26(χ²_min 1589.4, σ_eff 0.1189 dec) | results/e12_nd/summary.md |
| 실측 전달 시험 퍼짐 | 세트들이 예측한 Vth의 차이 3.6/1.2/1.3/8.2/3.6 mV(VD 0.1/0.2/0.5 V, L 2 µm VD 0.05/1 V; 각 세트의 전체 시뮬레이션 곡선) | results/e14/summary.md 2절, results/e13/summary.md |
| 실측 판별 예측 | t_ox 14 nm 133 mV, t_IWO 7 nm 376 mV, C–V 8.2 %. ND,eff 양 끝끼리만 47/126 mV | results/e13/summary.md |
| 보상 메커니즘 | 게이트 WF +170 mV 대 q·ΔN·t/C_ox +176 mV | results/e5_nd/summary.json의 θ로 계산 |
| DOS 조합 구분(E17) | 4%-W 조합 Δχ²_eff 30, Step-5 조합 13.5 | results/e14/summary.md 4절 |
| 시뮬레이터 차이 | Step-5 값이 측정보다 Vth +156/+208 mV | results/e10/summary.md |
| Step-5 두께 변형 Vth | t_IWO 7 nm −0.25 V, t_ox 14 nm −0.17 V(기본 −0.18 V), VD 1 V | results/crosscheck/step5_stack_vth.json |
| 결합 보정(E18), 쌍둥이 | VD 5개 160점, ND,eff 10^15.9–10^18.9에서 χ² 최대 0.644 | results/e18/summary.md |
| 결합 보정(E18), 실측 | VD 5개 93점, 10^16.85–10^19.05에서 Δχ²_eff 최대 2.25(χ²_min 1836.8, σ_eff 0.0935 dec). hold-out L 2 µm의 예측 Vth 차이 9.4/2.5 mV(VD 0.05/1 V) | results/e18/summary.md, results/e18_joint_meas/profile/summary.md |
| VD·L 시험의 해석적 근거 | 1차원 Pao–Sah: I(VG, VD) = qµ(W/L)∫_{VG−VD}^{VG} N_s(u) du, 2D Step-5와 Vth 1.7 mV 이내 | results/aux_paosah_1d/summary.md |
| 일반 규칙(E19), 점별 | σ_pred/σ 곡선 RMS: VD 0.02–3 V 0.25–0.41, L 0.5–10 µm 0.43–0.86, t_ox 5–20 nm 40.6–210, t_IWO 2–10 nm 92.2–487 | results/e19/summary.md |
| 일반 규칙(E19), 곡선 전체·비선형 | λ_max 수평 0.84–43 대 수직 1.0e5–1.9e7. E5 세트 10^15.9–10^17.9는 수평 시험 16개 χ² 합 0.07–0.11, 10^18.9는 시험별 최대 3.67(합 20.7; 보정만으로 맞춘 세트의 값이라 상한) | results/e19/summary.md |
| 재료 상수 ε_ch(E21) | ε_ch 15로 ND,eff 양 끝 재보정: χ² 1631.4/1608.4(기본 재료 1634.0/1606.3), 두 끝 Δχ²_eff 0.66, 기본 재료 χ²_min 기준 1.19/0.54 | results/e21/summary.md 1절, results/e14/summary.md 4절 |
| 시작점 의존성(E21) | 자유 보정 6개(v2·v3, 시작점 Step 5·4%-W·경계 중간값·E12 10^16.85) 모두 Δχ²_eff ≤ 1.10. N_TA 10^18.0–10^20.5(약 270배), 게이트 WF 182 mV, σ_GD 0.020–0.197 eV 차이. 전달 시험 Vth 예측 퍼짐 1.9–4.9 mV | results/e21/summary.md 2절 |
| 식별가능성 표(E20, 실측 결합 최적점) | 선형 1σ(σ_eff 0.0935 dec, Qf 고정). 보정 데이터만: ND,eff 7.3 dec(2곡선 8.0), N_TA 3.9 dec, 게이트 WF 1.21 eV. t_ox 14 nm·t_IWO 7 nm·C–V(0.5 %) 모두 추가: 이 점의 선형 1σ가 ND,eff 0.18 dec, 모든 descriptor 0.2 dec·15 meV 이내(C–V 2 %면 ND,eff 0.34, N_GD 0.45 dec). 점에 크게 좌우되는 국소값(2곡선 보정의 등가 점들에서 보정만 ND,eff 1σ 2.4–860 dec) | results/e20/summary.md, results/e21/summary.md |
| 공개 데모(E22) | φ_G–Q_f 등가 Q_f 4.5–10×10¹² cm⁻²(φ_G 4.20–4.55 eV), 등가 세트의 Vth 차이 VD·L 시험 0 mV 대 t_ox 14 nm 350 mV; 데모 컨트롤러 15회 Q_f 5.7–10×10¹²(t_ox 20 nm 오차 +5 ~ −499 mV) → t_ox 14 nm 곡선 1개 추가 시 5.62–5.67×10¹²; 계면 트랩 항 미지수면 N_tail 3.4 dec 등가, t_ch 분할만 가름(SS 차이 7 nm 33 mV/dec) | results/e22/summary.md, results/e22/summary.json |

수치를 글에 옮길 때는 반드시 위 파일에서 다시 확인합니다. 133/376 mV는 v2·v3 두 자유 보정 사이의 차이가 포함된 값이라 ND,eff 차이 때문이라고 쓰면 안 됩니다. 측정 대비 전달 시험 ΔVth는 디지타이즈 공백 때문에 근사값입니다(6장 "판독 신뢰도").

## 5. 다음 작업 구현 메모

### E18 결합 보정 profile — 완료(2026-10-24)

- 구현:
  - `e11_calibrate_measured.py`: `--vds`(보정 VD), `--max_gap`(측정 VG의 격자점만 사용). L1 hold-out에서는 보정에 쓴 VD를 뺍니다. J는 논문 정의대로 VD 0.05/1 V로만 계산합니다.
  - `e12_profile_measured.py`: `--vds`·`--max_gap`을 e11에 넘기고, 보고서 제목에 결합 VD를 적습니다.
  - `e5_profile_likelihood.py`: `--vds`(합성 데이터와 hold-out 목록에 반영). `e5b_merge_profiles.py`는 `vds`·`n_res`를 이어 받습니다.
  - `e18_report.py`: 판정표와 그림. L = 2 µm 예측 Vth의 세트 간 차이는 전체 시뮬레이션 곡선에서 계산합니다.
  - 저장된 결과로 보고서만 다시 만들면(README 10장) 기본값(VD 0.05/1 V) 결과가 바이트 단위로 같게 나옵니다. e11·e5를 새로 돌리면 summary.json에 `vds`(e11은 `max_gap`·`n_res_vd`도) 키가 추가되므로, 파일이 기존 것과 바이트 단위로 같지는 않습니다.
- 실행: profile 점마다 같은 ND,eff의 2곡선 결과에서 출발하는 독립 실행(README 10장 명령). 결과는 `results/e18_*`에만 둡니다. e14가 `results/e12_nd`·`e15_wta`·`e17_cross`를 읽어 χ²_min을 정하므로 섞이면 안 됩니다.
- 범위: 보정에 더한 것은 L = 1 µm의 VD 0.1/0.2/0.5 V 곡선뿐입니다. L = 2 µm 두 곡선은 hold-out으로 남겼습니다.
- 남은 보완(선택):
  - 실측 profile 점을 평가 40회로 다시 돌려 Δχ²_eff가 더 내려가는지 확인합니다(현재 기준점 대비 상한).
  - `fig_core_claim.py` (b)에 결합 보정 profile을 겹치는 개정은 아직 하지 않았습니다. 초록 그림(`fig_core_claim` 개정판, E18·E19 반영, README 6.4)과 함께 합니다.
  - L = 2 µm까지 넣은 7곡선 결합 보정: `CurveModel`은 기하 하나만 다루므로, L = 1·2 µm 두 모델의 잔차를 하나로 묶는 확장이 필요합니다. 실측 L = 2 µm 곡선도 공백이 있으므로 `--max_gap`을 함께 씁니다.

### E19 일반 규칙 — 완료(2026-10-24)

- 스크립트: `scripts/e19_design_rule.py`. 소자 목록은 `DEVICES`(라벨 → 기하, VD, 묶음)에 있습니다.
  - `--group lateral|vertical`: Step-5에서 조건별 감도 행(adjoint)을 계산합니다. 소자 하나가 끝날 때마다 `results/e19/rows_<group>.npz`와 `meta_<group>.json`에 저장하므로, 중간에 멈춰도 `--only <라벨>`로 이어 갈 수 있습니다.
  - `--nonlinear --values ... --tag ...`: E5 세트(보정 χ² < 3.84)를 수평 시험 조건에서 정방향으로 계산해 진짜와의 χ²를 냅니다(`nonlinear_<tag>.json`, 세트마다 저장).
  - `--report`: 시뮬레이션 없이 `summary.md`·`summary.json`·`e19_design_rule.png`를 만듭니다. 보정 행을 E3 행과, σ_pred/σ를 E14 값과 대조한 결과도 함께 적습니다.
- 보정 행(L 1 µm, VD 0.05/1 V)은 L1 소자 실행에서 함께 계산합니다(E3 `S_base`는 분석에 쓰지 않고 보고서에서 대조만 함).
- 실행 4개를 2개씩 병렬로 돌렸습니다(선형 45·33분, 비선형 78·58분). 비선형은 세트마다 수평 소자 다섯 개를 모두 다시 계산하므로 느립니다(L = 10 µm, 노드 24,897개).

### E20 식별가능성 표 — 완료(2026-10-24)

- 스크립트: `scripts/e20_identifiability_table.py`. E18 자유 보정(`results/e18_joint_meas/free`)의 θ에서 보정 잔차 행(e11과 같은 격자·바닥·`max_gap` 규칙)을 다시 만들고, t_ox 14 nm·t_IWO 7 nm(`sweep_with_grad`)·C–V(`cv_sweep`) 행을 계산합니다. 행은 `results/e20/rows.npz`, 보정 χ²와 C–V 점검은 `sim_meta.json`에 있습니다.
- `--report`는 시뮬레이션 없이 표를 다시 만들고, 쌍둥이 열은 `results/e3`·`e4`·`e7`의 행으로 계산합니다.
- 가정: 추천 I–V 시험에도 실측 보정의 σ_eff, C–V는 0.5 %(2 %) C_ox·A, Qf 고정(세 측정 모두에서는 Qf를 푼 열도 있음).

### E21 강건성 — 완료(2026-10-24)

- ε_ch = 15: `e11 --eps_ch 15 --fix log10_nd_eff=16.85`와 `19.05`를 돌렸습니다(같은 ND,eff의 E12 점에서 출발, 평가 25회, E16과 같은 방식). `e14`의 `robustness()`가 `results/e21_eps15/*`를 읽어 강건성 표에 넣습니다(두 끝 중 더 나은 쪽 기준).
- 다중 시작점: `results/e21_multistart/<step5|w4|mid|low_nd>`, 평가 40회.
  - `start_w4.json`(논문 4%-W 값)과 `start_mid.json`(경계 중간값)은 `theta_hat`만 담습니다. 그래서 Qf는 모든 실행에서 Step-5 값(1.59×10¹⁰ cm⁻²)입니다. `--start w4`를 쓰면 4%-W의 Qf(−1×10¹² cm⁻²)까지 바뀌어, Qf와 축퇴된 게이트 WF를 다른 실행과 비교할 수 없게 됩니다.
  - `low_nd`는 E12의 10^16.85 세트에서 ND,eff를 풀고 출발합니다.
- `scripts/e21_report.py`(시뮬레이션 없음): ε 표(같은 ND,eff의 기본 재료 E12 점과 χ²·descriptor 비교), 자유 보정 표(v2·v3 포함), 각 보정의 전체 시뮬레이션 곡선에서 구한 Vth 예측, 등가 도착점들의 퍼짐. χ²_min은 기본 재료의 모든 보정(E11 v2·v3, E12, E15, E17, 다중 시작점) 중 최저입니다.
- e12·e13·e14의 기준점(v3, χ²_min 1589.4)은 `results/e21_multistart`를 읽지 않습니다. 다중 시작점 보정이 v3보다 낮은 χ²에 도착하면 이 기준점들을 함께 바꿔야 합니다. 지금은 모두 v3보다 높습니다.

### 초록·논문 자료 (2026-10-24, 2026-10-25 갱신)

- `docs/paper/`: 핵심 수치표(`key_numbers.md`, `scripts/paper_numbers.py`로 저장된 결과에서 다시 계산, 상한은 올림), 참고문헌(`references.bib`, 확인하지 못한 서지는 `note`에 표시), 그림 목록과 서술 원칙.
- 초록 그림: `results/fig_kcs_abstract.png`(`scripts/fig_kcs_abstract.py`, 1행 3패널). v3(2026-10-25)은 (a) 워크플로 도식, (b) E19 설계 규칙, (c) E22 데모 끝점입니다. (b)는 E19, (c)는 E22의 `summary.json`만 읽습니다. 이전 판은 커밋 645f871에 있습니다.
- 초록 방향(v3): 기준 논문의 워크플로에 식별가능성 점검 단계(진단·재정의·설계·해소)를 더하는 것으로 씁니다. 논문이 틀렸다는 서술은 쓰지 않습니다.
- 초록·논문 원고는 저장소에 두지 않고 비공개로 관리합니다. 원고의 숫자는 `key_numbers.md`와 대조합니다.

### E22 공개 데모 모호성 분석 — 완료(2026-10-25)

- 스크립트: `scripts/e22_demo_ambiguity.py`(약 2분), 그림 `scripts/fig_e22_demo.py`. 결과는 `results/e22/summary.md`·`summary.json`·`fig_e22_demo.png`입니다.
- 환경: 데모를 설치해 import하므로 DEVSIM 환경과 별도인 Python 3.11 venv가 필요합니다(README 7장: `git clone https://github.com/gyulab/agentic-TCAD-calibration-public demo`, commit ccaa678, `requirements-lock.txt` + `pip install -e demo --no-deps`). 그림 스크립트는 어느 환경에서도 돌아갑니다(matplotlib만 사용).
- 데모 코드는 수정하지 않습니다. 두께 조건은 `GeoCondition`(데모 `Condition`의 하위 클래스)과 `GeoSimulator`(조건마다 `DeviceConstants`의 t_ox·t_ch를 바꾼 `CompactTFT`)로 넣었습니다. n_interface는 t_ox에 비례하게 바꿉니다(C_it 고정 가정). `free_cit=True`이면 `n_it`가 보정 파라미터가 됩니다.
- 데모의 `Workspace`·`Controller`·`RulePlanner`를 그대로 씁니다. 새로 만든 것은 `JointReachController`(활성화 판정에 joint-reach OR 조건)뿐이고, 기본 실행의 끝점이 비트 단위로 같은지 스크립트가 확인합니다(`identical_to_demo`).
- 결과는 결정적입니다. 같은 명령을 다시 돌리면 `summary.json`이 같습니다.
- 데모 저장소가 갱신되면 결과가 달라질 수 있으므로 commit ccaa678을 기준으로 둡니다.

### 문서 갱신

- 새 결과가 나오면 README 3장·6장, EXPERIMENT_LOG.md, `fig_core_claim.py`(E18·E19 반영 개정판)를 갱신합니다.
- 초록용 숫자는 README와 EXPERIMENT_LOG를 기준으로 합니다.

## 6. 알려진 함정

- **수렴:**
  - 긴 실행 중에 특정 세트가 바이어스 램프에서 실패할 수 있습니다. E1의 4%-W는 VG −1.475 V에서 실패해 `--only w4`로 따로 돌렸습니다.
  - E12의 10^19.2는 수렴하지 않아 10^19.05로 대신했습니다.
- **두꺼운 채널:** t_IWO 7 nm 소자는 Vth가 −0.77…−1.15 V(실측 세트, VD 1 V)입니다. 그래서 VG −3 V부터 계산합니다(`e13 --recompute "t_IWO 7 nm" --vgmin_iv -3.0`). 쌍둥이 Step-5는 −0.25 V입니다.
- **E13 두 번 실행:** 첫 실행은 모든 세트를 VG −1 V부터 계산하고, 두 번째 실행이 지정한 조건만 −3 V부터 다시 계산합니다.
- **χ²는 상한값:** 보정은 평가 횟수 한도(20–40회)에서 멈추는 경우가 대부분입니다. 그래서 χ²는 그 제약에서의 최솟값에 대한 상한입니다.
  - Δχ²_eff는 "지금까지의 모든 보정 중 최저 χ²"를 기준으로 하므로, 더 좋은 자유 보정이 생기면 모든 값이 바뀝니다. 예: v3가 생기면서 최대 0.77 → 최대 1.26. 기준점 v3 자체도 평가 한도(40회)와 경계(σ_GD = 0.02 eV)에서 멈춘 값입니다.
- **경로 의존성:** 자유 보정은 시작점과 반복 횟수에 따라 평평한 골짜기의 다른 점에 도착합니다(v2 대 v3, E21의 시작점 4개). 한 번의 보정 결과를 "그 데이터의 파라미터"로 쓰면 안 됩니다.
  - E21: 같은 시작점(Step 5)에서는 사실상 같은 점(v2)에 도착합니다. 시작점이 다르면 N_TA가 약 270배까지 다른 등가 점에 도착합니다(평가 40회 한도에서; 끝까지 수렴시키면 한 점으로 모이는지는 확인하지 않음). 낮은 ND,eff 세트에서 ND,eff를 풀고 출발하면 ND,eff가 10^16.85 → 10^18.29로 올라가지만 χ²는 σ_eff 기준 0.17만 줄어듭니다. 자유 보정 도착점의 ND,eff가 10^18.3–18.9에 모이는 것은 이 완만한 경사 때문으로 보이고, 식별을 뜻하지 않습니다.
  - 쌍둥이 E5의 10^19.4 점도 평가 20회와 30회 실행이 같은 χ²(187.79)에 다른 DOS 값을 줍니다. 저장된 값은 30회 실행(`e5_nd_up2`)입니다.
  - 평가 횟수 한도(`--max_nfev`)는 결과를 바꿉니다. README 10장 명령에 저장된 결과를 만든 값을 적었습니다(e11 기본 25, e12 기본 20, e5 기본 20).
- **경계에 걸린 파라미터:** 실측 profile의 10^16.85–10^17.85 세트는 E_GD = 0(하한), W_TA = 0.03 eV(하한)이고, v3와 E16 10^19.05는 σ_GD = 0.02 eV(하한), E17 4%-W 조합은 게이트 WF 3.70 eV(하한)입니다. 경계를 넓히면 각 세트의 χ²가 낮아질 수 있습니다. 기준점 χ²_min도 낮아질 수 있으므로 Δχ²_eff 자체는 상한이 아닙니다.
- **σ_eff:** 실측 보정에서 χ²/dof = 1이 되도록 키운 σ입니다. 실측 Δχ²는 항상 σ_eff로 환산해 3.84와 비교합니다.
- **DEVSIM 전역 상태:** 한 프로세스에서 TFT를 새로 build하면 이전 장치가 지워집니다(`reset_devsim`). 한 프로세스에서는 장치를 하나씩 다룹니다.
- **수렴 기준:** `absolute_error = 1e30`(상대오차만 사용), 아문턱 최대 반복 100회입니다. 민감도 선형계는 `ScaledLU`(equilibration + 반복 보정)를 씁니다. 일반 LU로 바꾸면 정확도가 5–6자리 떨어집니다.
- **E5 재실행:** README 10장의 네 명령을 그대로 쓰면 같은 결과가 나옵니다. 한 번에 기본값(평가 20회)으로 돌리면 구분 불가 영역은 마지막 자릿수, 10^19.15·10^19.4는 그보다 크게 달라질 수 있습니다.
- **논문 값 표기:** 2%-W Gaussian 중심은 Sec. III-A에서 Ec − 0.12 eV, Sec. IV에서 Ec − 0.122 eV입니다. 0.12 eV를 썼습니다.
- **판독 신뢰도:** VD 0.1/0.2/0.5 V 아문턱과 4%-W VD 0.05 V는 판독 점이 적습니다.
  - 디지타이즈한 전달 시험 곡선 다섯 개는 모두 100 nA 부근이 비어 있습니다. VD 0.1 V: −0.65 V 다음이 −0.1 V, VD 0.2 V: −0.65 V 다음이 0.45 V, VD 0.5 V: −0.85 V 다음이 0.35 V, L = 2 µm(VD 0.05·1 V): −0.55 V 다음이 −0.15 V. 이 곡선들의 "측정 Vth"는 보간값이고, VD 0.1/0.2/0.5 V는 사실상 정해지지 않습니다.
  - `holdout[k]["dVth_mV"]`(e11)는 모델을 측정 VG에서만 표본화해 같은 공백에 걸쳐 보간한 값입니다. 그래서 세트 간 비교에 쓰면 공백 때문에 값이 달라집니다(예: VD 0.5 V 36.8 mV 대 전체 곡선 1.3 mV). 세트 간 비교는 `holdout[k]["curve"]`와 `vg`로 전체 곡선의 Vth를 계산해서 합니다(e13 `full_curve_metrics`, e14 `vth_full_mV`, e18_report `full_vth`, `fig_core_claim.py` (c), `fig_measured.py` (d)).
  - L = 2 µm의 측정 대비 ΔVth는 표본화만으로 15–24 mV 달라집니다. 절대값을 쓸 때는 이 불확실성을 함께 적습니다.
  - 기본 E11 규칙은 측정 곡선을 격자에 선형 보간하므로, 이 곡선들을 보정에 넣으려면 `--max_gap 0.05`가 필요합니다(빈 구간 보간 금지).
- **E19 지표 읽기:**
  - σ_pred/σ는 점별 기준입니다. 곡선 전체에서 시험이 세트를 기각할 수 있는지는 λ_max(모든 방향)와 q(ND,eff profile 방향)로 봅니다. 수평 시험은 점별로는 1 미만이어도 곡선 전체로는 λ_max가 1을 넘을 수 있습니다(L = 2 µm 8.3–8.4). q는 나머지를 보정 곡선에만 다시 맞춘 값입니다. 시험까지 함께 다시 맞추면 profile Δχ²의 증가는 1/(ND,eff 1σ 비)²로 훨씬 작습니다.
  - χ²는 점 수에 비례합니다. 비선형 확인은 E5·E9와 같은 VG 0.1 V 간격이고, 선형 행은 0.05 V 간격입니다. 점 간격이 다르면 χ²를 직접 비교하지 않습니다.
  - 보정 곡선의 첫 점(10 pA 이상)은 VG −0.60 V입니다. L = 0.5 µm의 VD 1·3 V 곡선은 그보다 한 스텝 깊은 아문턱 점(−0.65 V)을 갖고, 그 점이 σ_pred/σ 최대값(3.6)을 만듭니다(보고서의 "max in window" 열). 비선형 확인의 0.1 V 격자에는 이 점이 없습니다.
- **E18 χ² 비교:** `--max_gap`을 쓰면 VD 0.05 V 잔차도 두 점(0.1, 0.2 V) 줄어듭니다. E18의 χ²·σ_eff(0.0935 dec)는 E11·E12(σ_eff 0.119 dec)와 직접 비교하지 않습니다.

## 7. 교차 확인 상태

- `results/crosscheck/`는 다른 시뮬레이터(예: 논문의 Sentaurus deck)에서 대안 세트를 돌려 보기 위한 묶음입니다.
  - `alternative_parameter_sets.csv`는 `scripts/make_crosscheck_package.py`로 다시 만듭니다.
  - `README.txt`에 정의, 두께·C–V 개념, 점검 절차가 있습니다.
- 확인 결과는 아직 없습니다(2026-10-24 기준).
- 결과가 오면 README 6.5의 판정 기준으로 반영합니다. 진행 경위 등 비공개 내용은 저장소 밖 노트에 있습니다.
