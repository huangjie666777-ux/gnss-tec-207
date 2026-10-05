# GNSS观测分析后端207

当前仓库：https://github.com/huangjie666777-ux/gnss-tec-207

纯后端服务：上传一份未压缩 RINEX 3.04 观测文件和一份 SP3-c 精密星历，
提供两条链路：

- `POST /position`：逐历元独立解算接收机 ECEF 位置与钟差（保留原能力）；
- `POST /ionosphere`：GPS 双频电离层监测，输出逐星斜向/垂直 TEC 与穿刺点。

## 环境

- Python 3.10.12 / FastAPI 0.115.12 / NumPy 2.2.6（全部在 .venv 中）

## 运行

    .venv/bin/python -m uvicorn main:app --host 127.0.0.1 --port 8000

## 使用

    curl -s -F "rinex=@examples/obs.rnx" -F "sp3=@examples/eph.sp3" \
         http://127.0.0.1:8000/position | python3 -m json.tool

返回每个历元：状态与失败原因、ECEF 米坐标、WGS84 经纬度与椭球高、
接收机钟差（秒）、使用/排除的卫星（含排除原因）、逐星残差（米）与 RMS。

## 电离层监测接口

	curl -s -X POST http://127.0.0.1:8000/ionosphere \
	  -F "rinex=@examples/obs.rnx" \
	  -F "sp3=@examples/eph.sp3" \
	  -F "lat_deg=30" -F "lon_deg=114" -F "height_m=50" \
	  -F "biases_ns=<examples/biases.json" | python3 -m json.tool

表单字段：`lat_deg/lon_deg/height_m` 为站点 WGS84 经纬高；
`biases_ns` 是 JSON 对象，逐星给出**合并码偏差（纳秒）**。
未出现在偏差表中的卫星每个历元标 `invalid`，不会静默消失。
返回逐星逐历元：GPS 时间、卫星号、弧 ID、状态
（`ok/short_arc/invalid/excluded`）、断弧或排除原因、斜向 TEC、
垂直 TEC、穿刺点地心经纬度与仰角。

## 生成合成示例数据

    .venv/bin/python examples/make_synthetic.py   # 写 examples/obs.rnx, examples/eph.sp3, examples/biases.json

## 测试

    .venv/bin/python -m pytest tests -q

## 处理规则与限制（共享）

- 仅 GPS 时间系统、GPS 卫星（G01–G32）、正常历元（flag 0）；
  接收机钟改正预应用（RCV CLOCK OFFS APPL != 0）、事件历元、非 3.04 版本、
  压缩文件、IONOSPHERIC CORR / TIME SYSTEM CORR / SYS SCALE FACTOR 等
  会改动原始信号的未支持修正记录均明确拒绝并指明文件与行号。
- 按头部观测类型顺序解析固定宽度字段，兼容 SYS / # / OBS TYPES 续行；
  空白或零伪距视为缺测。
- 必须同时声明 C1C、C2W、L1C、L2W；解析各观测量与相位 LLI 标志。
- 校验日期、卫星身份（G01–G32）、历元数量（上限 100）、非有限值与截断。
- SP3：位置 km→m，钟差 µs→s；零坐标与 999999.999999 缺失钟差被识别。
- 卫星位置用连续 8 节点 Lagrange 插值，钟差用相邻节点线性插值；
  不跨缺测节点、不外推，不可用的卫星给出排除原因。
- 每历元独立迭代最小二乘（位置 + 钟差），由伪距与卫星钟差推发射时刻，
  补偿信号飞行期间地球自转（Sagnac）；至少 4 颗有效卫星且几何满秩才求解；
  头部近似坐标只作初值，不沿用上一历元结果。
- 不足 4 星、秩亏或不收敛的历元返回原因，其余历元继续。

## 电离层处理规则

- **观测组合**：码差 P = C2W − C1C − c·b（b 为逐星合并码偏差，秒）；
  相位周转米 Lφ = L1C·λ1 − L2W·λ2，λ1、λ2 由 GPS L1/L2 频率与光速换算。
- **偏差符号**：输入偏差采用 TGD/DCB 类约定，即该偏差作为正延迟出现在
  观测码差中；因此从 P2−P1 **减去** c·b。请确认偏差产品符号约定与本服务一致。
- **分弧**：逐星分弧，任一观测量缺测、相位 LLI 非零、相邻历元间隔
  >120 s、或相邻相位差跳变 >1 m 即断弧；弧间**不共享**偏移。
  完整样本 <3 的弧输出 `short_arc` 并给出原因。
- **定级（弧内校准）**：至少 3 个完整样本的弧取
  b_arc = mean(P − Lφ)，再 STEC = (Lφ + b_arc)/K，
  K = 40.3e16·(1/f2² − 1/f1²)，单位 TECU；**负值保留**。
- **星历与几何**：复用 SP3 8 节点 Lagrange 插值在**接收历元**求卫星位置
  （补偿 Sagnac 地球自转），不跨缺失节点、不外推；否则该历元该星给出原因。
- **薄壳近似**：站星射线与半径 **6821 km** 地心球壳相交，站点径向为天顶；
  仰角 <10° 排除；垂直 TEC = 斜向 TEC × cos(z′)，z′ 为射线与壳面
  外法线夹角；穿刺点输出地心经纬度。薄壳不描述壳内电离层横向结构，
  低仰角映射误差较大。
- 缺测、无偏差、弧过短、仰角过低、星历不可用等记录都会返回并带原因。

## 精度边界

仅 C1C 伪距单点定位，无电离层/对流层/载波/多路径改正，
合成数据下亚米级；实际观测精度为米级，受大气延迟与观测噪声主导。

## 模块划分

- gnss_tec207/rinex.py — RINEX 3.04 解析与校验
- gnss_tec207/sp3.py — SP3-c 解析与单位换算
- gnss_tec207/interp.py — 位置 Lagrange / 钟差线性插值
- gnss_tec207/solver.py — 逐历元最小二乘定位
- gnss_tec207/arcs.py — 双频分弧、相位定级与斜向 TEC
- gnss_tec207/ionosphere.py — 薄壳穿刺几何与 TEC 流水线
- gnss_tec207/geodesy.py — WGS84 坐标转换与常数
- gnss_tec207/main.py — FastAPI 入口
