#!/bin/bash

# 1. 进入程序目录
cd /home/digitaldong/workplace/ProteinAnalyses/web_app

# 2. 检查程序是否已经在运行，如果在运行则杀掉旧进程
pkill -f "flask run"
pkill -f "lt --port 5000"

# 3. 后台启动 Flask 程序
echo "正在启动蛋白分析后端..."
export FLASK_APP=app.py
nohup python -m flask run --host=0.0.0.0 --port=5000 > flask.log 2>&1 &

# 4. 等待程序启动
sleep 3

# 5. 启动外网隧道 (使用随机二级域名)
echo "正在生成外网固定网址..."
nohup lt --port 5000 > lt.log 2>&1 &

# 6. 等待隧道建立并提取网址
sleep 5
URL=$(grep -o 'https://[a-zA-Z0-9.-]*\.loca\.lt' lt.log)

echo "------------------------------------------------"
echo "✅ 部署成功！"
echo "你的外网访问网址是: $URL"
echo "提示：打开网址时，如果看到 'Reminder' 页面，点击 'Click to Continue' 即可。"
echo "------------------------------------------------"
