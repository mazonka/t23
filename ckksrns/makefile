
MPIR=0
MPIR=1
MPIR=i


RT=bign$(MPIR)

include plat.mak


include mk_opts.mak

include plat_$(PLAT).mak

src1 = main.cpp mathut.cpp ckkselem.cpp ntt.cpp poly.cpp tmain.cpp \
	ckkshyb.cpp rns.cpp integer.cpp rns0.cpp

obj1 := $(src1:%.cpp=o/%.$(OEXT))

include mk_rules.mak

