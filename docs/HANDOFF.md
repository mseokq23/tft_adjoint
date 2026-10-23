# 인계 문서 — 새 세션에서 이어서 작업하기

이 문서는 작업을 새 환경이나 새 세션에서 정보 손실 없이 이어 가기 위한 것입니다. 연구 내용 해설은 [README.md](../README.md), 실험별 수치는 [EXPERIMENT_LOG.md](../EXPERIMENT_LOG.md)에 있습니다.

- 저장소: https://github.com/mseokq23/tft_adjoint (공개, 기본 브랜치 `main`)
- 기준 시점: 2026-10-23. 초록 마감: 2026-11-07(KCS 2027).
- 현재 상태: 핵심 가설 검증의 대부분이 끝났습니다(E0–E17). 남은 일은 E18–E21, 가치 판단(C1–C5), 초록 작성, 다른 시뮬레이터 교차 확인입니다(README 6장).

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
- 시뮬레이션 없이 다시 만들 수 있는 결과 목록은 README 10장에 있습니다. 모두 저장된 파일과 바이트 단위로 같게 재생성됨을 확인했습니다(2026-10-23).
- 지금까지의 실행 환경은 2 vCPU였고, 백그라운드 작업 하나가 2시간을 넘지 않게 나눠 돌렸습니다.
  - profile 계산은 점 단위로 이어 달리기(`e5 --init`, e12의 continuation)가 되므로 나눠 실행할 수 있습니다.
  - 동시에 두 프로세스까지 돌렸습니다.

## 4. 핵심 수치와 출처

| 주장 | 수치 | 파일 |
|---|---|---|
| Qf ≡ 게이트 WF | Fisher 고유값 비 2.2e-17, Qf 1e11 cm⁻² = WF 6.3 mV | results/e3/summary.md, results/e4/summary.md |
| donor ≡ acceptor + ND,eff | I–V ≤ 3.6e-15 dec, C–V ≤ 1.7e-14(두께 변형 ≤ 8.1e-15) | results/e8, results/e8_thickness |
| 쌍둥이 비식별 범위 | ND,eff 10^15.9–10^18.9에서 χ² ≤ 0.28, 95 % 상한 약 1.4e19 | results/e5_nd/summary.md |
| 쌍둥이 전달 시험 | max ΔVth ≤ 0.76 mV(10^19.4 세트도 13.7 mV) | results/e14/summary.md 1절 |
| 선형 기준 | 전달 시험 σ_pred/σ 0.25–0.43, 추천 시험 32–239 | results/e14/summary.md 3절 |
| 실측 비식별 범위 | ND,eff 10^16.85–10^19.05, Δχ²_eff ≤ 1.26(χ²_min 1589.4, σ_eff 0.1189 dec) | results/e12_nd/summary.md |
| 실측 전달 시험 퍼짐 | 3.2/4.6/5.5/3.6 mV(VD 0.5 V만 36.8 mV) | results/e14/summary.md 2절 |
| 실측 판별 예측 | t_ox 14 nm 133 mV, t_IWO 7 nm 376 mV, C–V 8.2 %. ND,eff 양 끝끼리만 47/126 mV | results/e13/summary.md |
| 보상 메커니즘 | 게이트 WF +170 mV 대 q·ΔN·t/C_ox +176 mV | results/e5_nd/summary.json의 θ로 계산 |
| DOS 조합 구분(E17) | 4%-W 조합 Δχ²_eff 30, Step-5 조합 13.5 | results/e14/summary.md 4절 |
| 시뮬레이터 차이 | Step-5 값이 측정보다 Vth +156/+208 mV | results/e10/summary.md |
| Step-5 두께 변형 Vth | t_IWO 7 nm −0.25 V, t_ox 14 nm −0.17 V(기본 −0.18 V), VD 1 V | results/crosscheck/step5_stack_vth.json |
| VD·L 시험의 해석적 근거 | 1차원 Pao–Sah: I(VG, VD) = qµ(W/L)∫_{VG−VD}^{VG} N_s(u) du, 2D Step-5와 Vth 1.7 mV 이내 | results/aux_paosah_1d/summary.md |

수치를 글에 옮길 때는 반드시 위 파일에서 다시 확인합니다. 133/376 mV는 v2·v3 두 자유 보정 사이의 차이가 포함된 값이라 ND,eff 차이 때문이라고 쓰면 안 됩니다.

## 5. 다음 작업 구현 메모

### E18 결합 보정 profile

- `scripts/e11_calibrate_measured.py`
  - 보정 VD `(0.05, 1.0)`이 다섯 곳에 직접 적혀 있습니다: 측정 선택 루프(102행), `CurveModel(... vds=...)`(112행), 보정 곡선 비교(140·142·147행).
  - `--vds`(쉼표 구분, 기본 `0.05,1.0`)를 추가해 이 다섯 곳에 쓰고, hold-out 목록(150행 `(0.1, 0.2, 0.5)`)에서는 보정에 쓴 VD를 뺍니다. 150행의 L = 2 µm 쌍 `(0.05, 1.0)`은 그대로 둡니다.
  - `summary.json`의 `cal` 키는 보정 VD별로 저장되므로, `e14_transfer_test_audit.py`의 `measured_sets()`(`("0.05", "1.0")`을 가정)도 함께 확인합니다.
- `scripts/e12_profile_measured.py`: `run_side()`가 e11에 `--max_nfev`만 넘깁니다. `--vds`를 받아 그대로 넘기게 합니다.
- `scripts/e5_profile_likelihood.py`(쌍둥이): 합성 데이터 생성과 `CurveModel(truth, names, vds=(0.05, 1.0), ...)`(약 86행)에 `--vds`를 반영합니다. `HOLDOUT`(37행)에서 보정에 쓴 VD를 뺍니다.
- **출력 폴더를 반드시 분리합니다**(`results/e18_joint_meas/`, `results/e18_joint_twin/`). 이유는 두 가지입니다.
  - `e14`의 `measured_sets()`는 `results/e12_nd/log10_nd_eff_*`를 모두 읽습니다.
  - `robustness()`는 `results/e15_wta/*`, `results/e17_cross/*`를 읽어 χ²_min을 정합니다.
  - 잔차 수가 다른(5 VD) 결합 보정의 χ²가 2곡선 분석의 χ²_min에 섞이면 안 됩니다.
- 순서:
  1. 결합 자유 보정(시작점: `results/e11_2W_v3/summary.json`, `--max_nfev 40`)
  2. 그 결과를 `--best`로 ND,eff profile(예: 18.35, 17.85, 16.85 / 19.05)
  3. 쌍둥이 2점
- 판정: 10배 이상 범위에서 Δχ²_eff < 3.84인지(σ_eff는 결합 보정의 값). VD 0.1/0.2/0.5 V 아문턱은 판독 점이 적다는 점을 함께 적습니다.

### E19 일반 규칙 (선형)

- 새 스크립트 `scripts/e19_design_rule.py`를 만듭니다.
- `e3_identifiability_oed_preview.py`의 `sweep_with_grad(params, geom, vds, vgs)`와 `log_rows(vg, i, G)`로 Step-5에서 조건별 감도 행을 계산합니다.
- `e14`의 `linear_predictive()`와 같은 방식으로 σ_pred/σ를 냅니다. 보정 행은 `results/e3/rows.npz`의 `S_base`, Qf 열은 뺍니다.
- 후보:
  - 수평·바이어스: `Geometry(L_um=0.5/5/10)`, VD 0.02–3 V
  - 수직 구조: `Geometry(t_gox_nm=5–20)`, `Geometry(t_iwo_nm=2–10)`
- 출력: `results/e19/`(rows.npz, summary.md, 그림)
- 주의: VD 3 V와 L 0.5 µm는 바이어스 램프 수렴을 먼저 확인합니다.
- 해석적 기준: `scripts/aux_paosah_1d.py`(1차원 Pao–Sah)에서는 VD·L 시험이 정확히 정보를 더하지 않습니다. E19 결과를 이 기준과 비교하면, 2D 모델에서 벗어나는 부분(접촉·겹침 영역, 단채널 효과)이 드러납니다.

### E20 식별가능성 표

- E18 결합 최적점의 θ로 E3(I–V), E4·E7(두께, C–V) 방식의 감도 행을 다시 계산합니다.
- Fisher 1σ(σ_eff)를 "보정만 / + t_ox 14 nm / + t_IWO 7 nm / + C–V / 모두"로 표 하나에 정리합니다.

### E21 강건성

- ε_ch = 15: `e11 --eps_ch 15 --fix log10_nd_eff=16.85`와 `19.05`를 돌립니다(E16과 같은 방식). `e14`의 `robustness()`에 `results/e21_eps15/*` 행을 추가합니다.
- 다중 시작점: `e11 --start`에 서로 다른 시작점 4개(예: `step5`, `w4`, 일반 초기값 JSON, v3)를 줘서 자유 보정 도착점을 비교합니다.

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
  - Δχ²_eff는 "지금까지의 모든 보정 중 최저 χ²"를 기준으로 하므로, 더 좋은 자유 보정이 생기면 모든 값이 바뀝니다. 예: v3가 생기면서 ≤ 0.77 → ≤ 1.26.
- **경로 의존성:** 자유 보정은 시작점과 반복 횟수에 따라 평평한 골짜기의 다른 점에 도착합니다(v2 대 v3). 한 번의 보정 결과를 "그 데이터의 파라미터"로 쓰면 안 됩니다.
  - 쌍둥이 E5의 10^19.4 점도 평가 20회와 30회 실행이 같은 χ²(187.79)에 다른 DOS 값을 줍니다. 저장된 값은 30회 실행(`e5_nd_up2`)입니다.
  - 평가 횟수 한도(`--max_nfev`)는 결과를 바꿉니다. README 10장 명령에 저장된 결과를 만든 값을 적었습니다(e11 기본 25, e12 기본 20, e5 기본 20).
- **경계에 걸린 파라미터:** 실측 profile의 10^16.85–10^17.85 세트는 E_GD = 0(하한), W_TA = 0.03 eV(하한)이고, v3와 E16 10^19.05는 σ_GD = 0.02 eV(하한), E17 4%-W 조합은 게이트 WF 3.70 eV(하한)입니다. 경계를 넓히면 χ²가 낮아질 수 있습니다(Δχ²_eff는 상한).
- **σ_eff:** 실측 보정에서 χ²/dof = 1이 되도록 키운 σ입니다. 실측 Δχ²는 항상 σ_eff로 환산해 3.84와 비교합니다.
- **DEVSIM 전역 상태:** 한 프로세스에서 TFT를 새로 build하면 이전 장치가 지워집니다(`reset_devsim`). 한 프로세스에서는 장치를 하나씩 다룹니다.
- **수렴 기준:** `absolute_error = 1e30`(상대오차만 사용), 아문턱 최대 반복 100회입니다. 민감도 선형계는 `ScaledLU`(equilibration + 반복 보정)를 씁니다. 일반 LU로 바꾸면 정확도가 5–6자리 떨어집니다.
- **E5 재실행:** 저장된 E5는 네 번에 나눠 계산해 합친 것이라 다시 돌리면 마지막 자릿수가 다를 수 있습니다.
- **논문 값 표기:** 2%-W Gaussian 중심은 Sec. III-A에서 Ec − 0.12 eV, Sec. IV에서 Ec − 0.122 eV입니다. 0.12 eV를 썼습니다.
- **판독 신뢰도:** VD 0.1/0.2/0.5 V 아문턱과 4%-W VD 0.05 V는 판독 점이 적습니다. 실측 VD 0.5 V의 세트 간 퍼짐(36.8 mV)이 큰 것도 이 때문일 수 있습니다.

## 7. 교차 확인 상태

- `results/crosscheck/`는 다른 시뮬레이터(예: 논문의 Sentaurus deck)에서 대안 세트를 돌려 보기 위한 묶음입니다.
  - `alternative_parameter_sets.csv`는 `scripts/make_crosscheck_package.py`로 다시 만듭니다.
  - `README.txt`에 정의, 두께·C–V 개념, 점검 절차가 있습니다.
- 확인 결과는 아직 없습니다(2026-10-23 기준).
- 결과가 오면 README 6.5의 판정 기준으로 반영합니다. 진행 경위 등 비공개 내용은 저장소 밖 노트에 있습니다.
