# 초록·논문용 자료

학회 초록과 이후 전체 논문을 쓸 때 쓰는 자료를 모았습니다. 원고 자체는 저장소 밖에서 관리하고, 여기에는 다시 만들 수 있는 그림·숫자표·참고문헌만 둡니다.

## 파일

| 파일 | 내용 | 만드는 방법 |
|---|---|---|
| `key_numbers.md`, `key_numbers.json` | 초록·논문에 쓰는 핵심 수치와 각 수치의 출처 파일 | `python scripts/paper_numbers.py` (시뮬레이션 없음, 저장된 결과에서 다시 계산) |
| `references.bib` | 참고문헌. 쪽수·DOI를 확인하지 못한 항목은 `note`에 적어 둠 | 수동 관리 |
| `../../results/fig_kcs_abstract.png` | 1쪽 초록용 1행 3패널 그림 | `python scripts/fig_kcs_abstract.py` |

결과 파일이 바뀌면 `paper_numbers.py`와 그림 스크립트를 다시 돌리고, 원고의 숫자를 `key_numbers.md`와 대조합니다.

## 초록 그림 (`results/fig_kcs_abstract.png`)

- (a) ND,eff profile. 쌍둥이(σ = 0.02 dec)와 디지타이즈 실측(σ_eff), 각각 보정 곡선 2개(VD 0.05/1 V)와 드레인 바이어스 곡선 5개(E5, E18, E12). 95 % 기준과 논문의 2%-W·4%-W ND,eff 값을 함께 표시합니다.
- (b) 실측에 똑같이 맞는 7세트(E12 profile 점과 자유 보정 v2·v3)의 예측 Vth. 논문 전달 시험(각 세트의 전체 시뮬레이션 곡선)과 t_ox 14 nm·t_IWO 7 nm 소자(E13, VD 1 V)를 시험마다 중앙값 세트 기준으로 그렸습니다. 가장 큰 두께 시험 차이는 DOS가 다른 v2·v3 사이에서 나옵니다(ND,eff 양 끝만 비교하면 47/126 mV).
- (c) 시험 32개의 예측 σ ÷ 시험 잡음(E19, 쌍둥이, 선형, 곡선 RMS). 드레인 바이어스(VD 0.02–3 V)·채널 길이(L 0.5–10 µm) 시험은 1 미만이고, 게이트 절연막(5–20 nm)·채널(2–10 nm) 두께 시험은 41–487입니다.

## 논문에 쓸 수 있는 다른 그림

| 그림 | 내용 |
|---|---|
| `results/fig_core_claim.png` | 핵심 주장 2×2: 같은 실측 곡선, profile, 예측 Vth 퍼짐, σ/잡음 |
| `results/fig_measured.png` | 실측 비교 2×2 |
| `results/fig_abstract.png` | 쌍둥이 2×2 |
| `results/e18/e18_profile.png` | 결합 보정(드레인 바이어스 곡선 5개) profile |
| `results/e19/e19_design_rule.png` | E19 일반 규칙 6패널 |
| `results/e13/e13_predictions.png` | 실측 등가 세트의 두께·C–V 예측 |
| `results/aux_paosah_1d/` | 1차원 Pao–Sah 해석(VD·L 시험이 정보를 더하지 않는 이유) |
| `results/e22/fig_e22_demo.png` | 기준 논문 저자의 공개 데모에서의 확인(E22): Q_f profile, 데모 컨트롤러 15회 끝점, 계면 트랩 항을 미지수로 둘 때의 N_tail profile |

## 서술 원칙

README 6.4와 5장을 따릅니다.

- "최초(first)", "exact gradient", "LLM 에이전트보다 우수" 같은 표현은 쓰지 않습니다.
- 비유일성 자체, 두께 분할·C–V를 이용한 분리, donor/acceptor 항등은 알려진 것으로 인용합니다. 기여는 전달 시험 논리의 정량적 반증과 수평/수직 시험 설계 규칙입니다.
- 범위(같은 descriptor 클래스, 이 오픈소스 구현, 디지타이즈 데이터)를 밝힙니다. 실측 결론은 모델 오차를 포함한 σ_eff 기준의 상대 비교입니다.
- E20의 1σ는 선형화한 점에 크게 좌우되는 국소값입니다. 절대값보다 측정을 더했을 때 줄어드는 폭으로 씁니다.
- 133/376 mV는 v2·v3 차이가 포함된 값입니다. ND,eff 양 끝 차이(47/126 mV)와 구분합니다.
