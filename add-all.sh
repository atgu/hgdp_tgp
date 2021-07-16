for POLICY in "max-10" "max-50" "max-100" "max-500"; do
    gcloud dataproc --project neurogap-analysis autoscaling-policies import ${POLICY} --source \
      <(curl https://raw.githubusercontent.com/hail-is/dataproc-autoscaling-configs/master/${POLICY}-preemptibles.yaml)
done
