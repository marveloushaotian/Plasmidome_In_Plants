#!/bin/bash

# 自动检测所有chunk文件并并行运行
max_parallel=10  # 同时运行的任务数，可以根据资源调整

for chunk_file in eggnog_input/eggnog_inputs_prokka.part_*; do
    # 控制并发数
    while [ $(jobs -r | wc -l) -ge $max_parallel ]; do
        sleep 5  # 等待5秒再检查
    done
    
    echo "启动处理: $chunk_file"
    emapper.py -m diamond --no_annot --no_file_comments --cpu 20 \
        -i "$chunk_file" -o "$chunk_file" \
        --data_dir ../../eggnog_db/ -d none --itype proteins --override \
        1>"${chunk_file}.log" 2>&1 &
done

# 等待所有任务完成
wait
echo "所有diamond搜索完成！"
