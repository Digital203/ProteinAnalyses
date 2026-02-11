# 蛋白分析（通用单蛋白模板）

该目录提供一个可复用的单蛋白 PDB 分析脚本，支持快速替换任意蛋白进行统一科研可视化输出。

## 依赖

```bash
pip install numpy matplotlib
```

## 用法

在 `蛋白分析` 目录运行：

```bash
python protein_analysis.py \
  --pdb ../52.pdb \
  --chain A \
  --output-dir ./results_52 \
  --prefix 52
```

## 输出文件

- `*_summary.json`：结构与统计摘要（长度、序列、氨基酸组成、回转半径、分数统计）
- `*_3d_backbone.png`：CA 骨架三维图（按分数着色）
- `*_residue_score.png`：每残基分数曲线
- `*_ca_distance_map.png`：CA 距离矩阵图
- `*_ramachandran.png`：Ramachandran 图
- `*_overview.png`：总览图（分数曲线 + 氨基酸组成）

## 通用接口说明

- 更换蛋白：只改 `--pdb` 路径
- 指定链：改 `--chain`（不指定则分析所有链）
- 修改输出目录/前缀：改 `--output-dir` 和 `--prefix`
- 修改图清晰度：改 `--dpi`
