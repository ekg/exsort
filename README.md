# exsort

## tests

```
cat /dev/urandom | head -c $(echo '8 * 10000000' | bc) >y 
cp y z && time ./exsort -r 8 -k 8 -t 4 z           
```

