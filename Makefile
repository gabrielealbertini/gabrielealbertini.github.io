.phony: cv all

all: cv

cv:
	cp ../curriculum_vitae/gabriele_albertini_vitae.pdf .
	git add gabriele_albertini_vitae.pdf
	git commit -m "update cv"
	git push
