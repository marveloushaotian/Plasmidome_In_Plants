#!/bin/bash
# Safe pipeline runner with proper Java memory settings

# Set Nextflow Java heap memory
export NXF_OPTS='-Xms2g -Xmx8g'

# Set number of worker threads to avoid overwhelming the system
export NXF_EXECUTOR_PERJOBTMP=true

echo "=========================================="
echo "Pipeline Runner"
echo "=========================================="
echo "Java heap: 2GB-8GB"
echo "Nextflow version: $(nextflow -version 2>&1 | head -1)"
echo "=========================================="
echo ""

# Run the pipeline
nextflow run main.nf \
  --run_gtdbtk true \
  -profile standard \
  -resume \
  -with-dag flowchart.svg \
  -with-timeline timeline.html \
  -with-report report.html \
  -with-trace trace.txt

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Pipeline completed successfully!"
else
    echo "✗ Pipeline failed with exit code: $EXIT_CODE"
    echo ""
    echo "Check logs:"
    echo "  tail -100 .nextflow.log"
fi
echo "=========================================="

exit $EXIT_CODE

