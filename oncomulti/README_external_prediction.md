# 外部预测集功能使用说明

## 功能概述

本系统现在支持将指定的样本作为外部预测集，用于模型训练完成后的外部验证。这样可以更好地评估模型的泛化能力。

## 主要特性

1. **自动数据分割**: 将指定的10个样本从训练数据中分离出来作为外部预测集
2. **完整的外部验证**: 包括性能指标计算、可视化图表生成
3. **详细预测结果**: 保存每个样本的预测结果和概率
4. **独立的外部预测脚本**: 可以单独运行外部预测

## 外部预测集样本

默认的外部预测集包含以下样本（共9个可用）：

- OSCC115
- OSCC11
- OSCC62
- OSCC29
- OSCC125
- OSCC120
- OSCC33
- OSCC130
- OSCC88


## 使用方法

### 1. 完整训练和外部验证

运行主程序，系统会自动进行数据分割、模型训练和外部验证：

```bash
python main.py
```

这将生成：

- 训练集上的模型性能指标
- 外部预测集上的性能指标
- 各种可视化图表（ROC曲线、混淆矩阵等）
- 详细的预测结果文件

### 2. 仅进行外部预测

如果模型已经训练完成，可以单独运行外部预测：

```bash
python external_prediction.py
```

### 3. 自定义外部预测集

可以修改 `external_prediction.py` 中的 `external_samples` 参数来指定不同的外部预测集样本。

## 输出文件

### 训练阶段输出

- `results/model_report.json`: 训练集评估报告
- `results/external_prediction_report.json`: 外部预测集评估报告
- `results/external_prediction_data.csv`: 外部预测集详细数据
- `plots/roc_curve.png`: 训练集ROC曲线
- `plots/external_roc_curve.png`: 外部预测集ROC曲线
- `plots/confusion_matrix.png`: 训练集混淆矩阵
- `plots/external_confusion_matrix.png`: 外部预测集混淆矩阵

### 外部预测阶段输出

- `results/external_prediction_results.csv`: 详细预测结果
- `plots/external_prediction_roc.png`: 外部预测ROC曲线
- `plots/external_prediction_confusion_matrix.png`: 外部预测混淆矩阵

## 数据分割详情

- **训练集**: 222个样本（健康对照=129, 口腔鳞癌=93）
- **外部预测集**: 9个样本（健康对照=0, 口腔鳞癌=9）
- **总样本数**: 231个样本
- **数据完整性**: 训练集和外部预测集完全分离，无数据泄露

## 性能指标

系统会计算以下性能指标：

- F1分数（主要指标）
- AUC
- 准确率
- 精确率
- 召回率（敏感性）
- 特异性

## 注意事项

1. 外部预测集样本必须在原始数据中存在
2. 系统会自动检查样本是否存在，并报告缺失的样本
3. 外部预测集主要用于评估模型泛化能力，不参与模型训练
4. 建议外部预测集包含不同类别的样本以获得更全面的评估
