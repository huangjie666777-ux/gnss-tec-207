# GNSS观测分析后端207

当前仓库：https://github.com/huangjie666777-ux/gnss-tec-207

纯后端服务：上传一份未压缩 RINEX 3.04 观测文件和一份 SP3-c 精密星历，
逐历元独立解算接收机 ECEF 位置与钟差。

## 环境

- Python 3.10.12 / FastAPI 0.115.12 / NumPy 2.2.6（全部在 .venv 中）

## 运行

    .venv/bin/python -m uvicorn main:app --host 127.0.0.1 --port 8000

## 使用

    curl -s -F "rinex=@examples/obs.rnx" -F "sp3=@examples/eph.sp3" \
         http://127.0.0.1:8000/position | python3 -m json.tool

返回每个历元：状态与失败原因、ECEF 米坐标、WGS84 经纬度与椭球高、
接收机钟差（秒）、使用/排除的卫星（含排除原因）、逐星残差（米）与 RMS。

## 生成合成示例数据

    .venv/bin/python examples/make_synthetic.py   # 写 examples/obs.rnx, examples/eph.sp3

## 测试

    .venv/bin/python -m pytest tests -q

## 处理规则与限制

- 仅 GPS 时间系统、GPS 卫星、C1C 伪距、正常历元（flag 0）；
  接收机钟改正预应用（RCV CLOCK OFFS APPL != 0）、事件历元、非 3.04 版本、
  压缩文件等明确拒绝并指明文件与行号。
- 按头部观测类型顺序解析固定宽度字段，兼容 SYS / # / OBS TYPES 续行；
  空白或零伪距视为缺测。
- 校验日期、卫星身份、历元数量（上限 100）、非有限值与截断。
- SP3：位置 km→m，钟差 µs→s；零坐标与 999999.999999 缺失钟差被识别。
- 卫星位置用连续 8 节点 Lagrange 插值，钟差用相邻节点线性插值；
  不跨缺测节点、不外推，不可用的卫星给出排除原因。
- 每历元独立迭代最小二乘（位置 + 钟差），由伪距与卫星钟差推发射时刻，
  补偿信号飞行期间地球自转（Sagnac）；至少 4 颗有效卫星且几何满秩才求解；
  头部近似坐标只作初值，不沿用上一历元结果。
- 不足 4 星、秩亏或不收敛的历元返回原因，其余历元继续。

## 精度边界

仅 C1C 伪距单点定位，无电离层/对流层/载波/多路径改正，
合成数据下亚米级；实际观测精度为米级，受大气延迟与观测噪声主导。

## 模块划分

- gnss_tec207/rinex.py — RINEX 3.04 解析与校验
- gnss_tec207/sp3.py — SP3-c 解析与单位换算
- gnss_tec207/interp.py — 位置 Lagrange / 钟差线性插值
- gnss_tec207/solver.py — 逐历元最小二乘定位
- gnss_tec207/geodesy.py — WGS84 坐标转换与常数
- gnss_tec207/main.py — FastAPI 入口
