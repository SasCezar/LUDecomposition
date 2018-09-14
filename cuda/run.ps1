$blocks = 2, 4, 6, 8
$threads = 256, 512, 1024
For ($i=1; $i -le 10 ; $i++) {
	For ($s=1000; $s -le 9000; $s = $s + 1000) {
		Foreach ($b in $blocks) {
			Foreach ($t in $threads) {
				./x64/Release/cuda.exe $s $t $b
			}
		}
	}
}
