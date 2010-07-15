#!/bin/bash

postfix=-impi

cd $PWD

# comment in everything 
sed -i 's/^PARALLEL += $(FDEF)-DISEND_IRECV/#PARALLEL += $(FDEF)-DISEND_IRECV/' Makefile
sed -i 's/^PARALLEL += $(FDEF)-DSENDRECV/#PARALLEL += $(FDEF)-DSENDRECV/' Makefile
sed -i 's/^PARALLEL += $(FDEF)-DMPI_SUBARRAY/#PARALLEL += $(FDEF)-DMPI_SUBARRAY/' Makefile
sed -i 's/^PARALLEL += $(FDEF)-DUSE_ADCL/#PARALLEL += $(FDEF)-DUSE_ADCL/' Makefile

sed -i 's/^#PARALLEL += $(FDEF)-DVERBOSE/PARALLEL += $(FDEF)-DVERBOSE/' Makefile

# build isir
echo "... building isir"
sed -i 's/#PARALLEL += $(FDEF)-DISEND_IRECV/PARALLEL += $(FDEF)-DISEND_IRECV/' Makefile
make clean && make
mv lbc lbc-no-adcl-isir$postfix
sed -i 's/PARALLEL += $(FDEF)-DISEND_IRECV/#PARALLEL += $(FDEF)-DISEND_IRECV/' Makefile

# build sr 
echo "... building sr"
sed -i 's/#PARALLEL += $(FDEF)-DSENDRECV/PARALLEL += $(FDEF)-DSENDRECV/' Makefile
make clean && make
mv lbc lbc-no-adcl-sr$postfix
sed -i 's/PARALLEL += $(FDEF)-DSENDRECV/#PARALLEL += $(FDEF)-DSENDRECV/' Makefile

# build subarr
echo "... building subarr"
sed -i 's/#PARALLEL += $(FDEF)-DMPI_SUBARRAY/PARALLEL += $(FDEF)-DMPI_SUBARRAY/' Makefile
make clean && make
mv lbc lbc-no-adcl-subarr$postfix
sed -i 's/PARALLEL += $(FDEF)-DMPI_SUBARRAY/#PARALLEL += $(FDEF)-DMPI_SUBARRAY/' Makefile

# build adcl
echo "... building adcl"
sed -i 's/#PARALLEL += $(FDEF)-DUSE_ADCL/PARALLEL += $(FDEF)-DUSE_ADCL/' Makefile
make clean && make
mv lbc lbc-with-adcl$postfix
sed -i 's/PARALLEL += $(FDEF)-DUSE_ADCL/#PARALLEL += $(FDEF)-DUSE_ADCL/' Makefile

echo "done"

