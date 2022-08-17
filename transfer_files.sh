#!/bin/bash
for i in {0 ... 1999}
	do sshpass -p "clkhin13" scp -r d$i khalil@192.168.1.23:/home/khalil/Repos/mvm_paper/data/strut_fea/opt_manual/d$i
done
