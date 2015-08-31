genome=../ath/TAIR10_Chr.all.fasta
mkfifo r1
mkfifo r2
pairs join r1 r2 >ath_fake_novar_1M.ilfq &
wgsim -S 3301 -r 0.0 -N 1000000 -1 101 -2 101 -d 300 -s 100 $genome r1 r2

pairs join r1 r2 >ath_fake_1e-3var_1M.ilfq &
wgsim -S 3301 -r 0.001 -N 1000000 -1 101 -2 101 -d 300 -s 100 $genome r1 r2 >ath_fake_1e-3var_1M.var

pairs join r1 r2 >ath_fake_1e-2var_1M.ilfq &
wgsim -S 3301 -r 0.01 -N 1000000 -1 101 -2 101 -d 300 -s 100 $genome r1 r2 >ath_fake_1e-2var_1M.var
rm -v r1 r2
