#!/bin/bash

# 带并发控制的批量处理eggNOG-mapper注释脚本

# =============配置参数=============
INPUT_DIR="Intermediate_files_diamond"
CPU_PER_TASK=10          # 每个任务使用的CPU数
TOTAL_CPU=120            # 服务器总CPU数
MAX_CONCURRENT=$((TOTAL_CPU / CPU_PER_TASK))  # 最大并发任务数

echo "=== 批量eggNOG-mapper注释任务 ==="
echo "服务器总CPU: $TOTAL_CPU"
echo "每任务CPU: $CPU_PER_TASK" 
echo "最大并发数: $MAX_CONCURRENT"
echo "输入目录: $INPUT_DIR"
echo "=================================="

# 检查输入文件夹是否存在
if [ ! -d "$INPUT_DIR" ]; then
    echo "❌ 错误: 文件夹 $INPUT_DIR 不存在!"
    exit 1
fi

# 创建文件列表
FILES=()
for file in "$INPUT_DIR"/*.emapper.seed_orthologs; do
    if [ -f "$file" ]; then
        FILES+=("$file")
    fi
done

if [ ${#FILES[@]} -eq 0 ]; then
    echo "❌ 没有找到 .emapper.seed_orthologs 文件"
    exit 1
fi

echo "📋 找到 ${#FILES[@]} 个文件需要处理"
echo ""

# 存储正在运行的任务信息
declare -A running_tasks  # 关联数组：PID -> 文件名
current_file_index=0
completed_count=0
total_files=${#FILES[@]}

# 函数：启动单个任务
start_task() {
    local file="$1"
    local filename=$(basename "$file")
    local basename_no_ext=$(basename "$file" .emapper.seed_orthologs)
    local output_name="${basename_no_ext}_eggnog"
    local log_file="annotation_${basename_no_ext}.log"
    
    echo "🚀 启动任务: $filename"
    echo "   输出前缀: $output_name"
    echo "   日志文件: $log_file"
    
    # 启动任务
    emapper.py --annotate_hits_table "$file" \
               --no_file_comments \
               -o "$output_name" \
               --cpu $CPU_PER_TASK \
               --data_dir /dev/shm/eggnog_db \
               1>"$log_file" 2>&1 &
    
    local pid=$!
    running_tasks[$pid]="$filename"
    
    echo "   ✅ 已启动 PID: $pid"
    echo "   📊 运行中任务: ${#running_tasks[@]}/$MAX_CONCURRENT"
    echo ""
}

# 函数：检查并清理已完成的任务
check_completed_tasks() {
    local completed_pids=()
    
    # 检查每个运行中的任务
    for pid in "${!running_tasks[@]}"; do
        if ! kill -0 "$pid" 2>/dev/null; then
            # 任务已完成
            completed_pids+=($pid)
            completed_count=$((completed_count + 1))
            echo "✅ 任务完成: ${running_tasks[$pid]} (PID: $pid)"
            echo "   📈 进度: $completed_count/$total_files"
        fi
    done
    
    # 从running_tasks中移除已完成的任务
    for pid in "${completed_pids[@]}"; do
        unset running_tasks[$pid]
    done
    
    echo "   📊 当前运行: ${#running_tasks[@]}/$MAX_CONCURRENT"
    echo ""
}

# 函数：等待直到有可用的槽位
wait_for_available_slot() {
    while [ ${#running_tasks[@]} -ge $MAX_CONCURRENT ]; do
        echo "⏳ 已达到最大并发数 ($MAX_CONCURRENT)，等待任务完成..."
        sleep 10  # 每10秒检查一次
        check_completed_tasks
    done
}

# 主循环：处理所有文件
echo "🔄 开始批量处理..."
echo ""

while [ $current_file_index -lt $total_files ]; do
    # 检查并清理已完成的任务
    check_completed_tasks
    
    # 等待可用槽位
    wait_for_available_slot
    
    # 启动新任务
    file="${FILES[$current_file_index]}"
    start_task "$file"
    
    current_file_index=$((current_file_index + 1))
    
    # 显示总体进度
    echo "📋 总体进度: 已启动 $current_file_index/$total_files, 已完成 $completed_count/$total_files"
    echo "----------------------------------------"
done

# 等待所有剩余任务完成
echo "🏁 所有任务已启动，等待剩余任务完成..."
echo ""

while [ ${#running_tasks[@]} -gt 0 ]; do
    echo "⏳ 等待剩余 ${#running_tasks[@]} 个任务完成..."
    echo "   运行中的任务:"
    for pid in "${!running_tasks[@]}"; do
        echo "     - PID $pid: ${running_tasks[$pid]}"
    done
    echo ""
    
    sleep 15  # 每15秒检查一次
    check_completed_tasks
done

# 最终结果
echo ""
echo "🎉 =============== 任务完成 ==============="
echo "📊 总处理文件: $total_files"
echo "✅ 成功完成: $completed_count"
echo ""

# 检查生成的文件
echo "📁 检查生成的结果文件:"
result_count=$(ls *_eggnog.emapper.annotations 2>/dev/null | wc -l)
if [ $result_count -gt 0 ]; then
    echo "✅ 找到 $result_count 个注释结果文件:"
    ls -1 *_eggnog.emapper.annotations 2>/dev/null | head -10
    if [ $result_count -gt 10 ]; then
        echo "   ... 还有 $((result_count - 10)) 个文件"
    fi
else
    echo "❌ 未找到注释结果文件，请检查日志"
fi

echo ""
echo "📋 查看详细结果:"
echo "   - 所有结果文件: ls -la *_eggnog.*"
echo "   - 查看日志错误: grep -i error annotation_*.log"
echo "   - 查看任务完成状态: grep -i 'done\|finished\|completed' annotation_*.log"
echo "==============================================="
