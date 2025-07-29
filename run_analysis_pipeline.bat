@echo off

:: 激活虚拟环境
call venv\Scripts\activate.bat

:: 生成适体库
python generate_aptamer_libraries.py

:: 预测二级结构
python predict_secondary_structures.py aptamer_libraries/all_aptamers_combined.fasta -o secondary_structures.txt -t 8

:: 分析结构特征
python analyze_aptamer_structures.py

:: 退出虚拟环境
deactivate
echo 完整分析流程已完成!
endlocal
