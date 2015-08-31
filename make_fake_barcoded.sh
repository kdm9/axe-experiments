
python fakebcd.py ath_fake_novar_1M.ilfq >ath_fake_novar_1M_barcoded.ilfq 2>ath_fake_novar_1M_barcoded.barcodes

paste <(cut -f 1 ath_fake_novar_1M_barcoded.barcodes |grep -v None) <(for i in {1..10}; do echo $i; done) > fake_barcodes.axe
