find . -type f -name '*:*' | while read -r file; do
  mv -- "$file" "${file//:/-}"
done