all:
	python weather.py

fresh: clean
	python weather.py

cleancache:
	rm -rf __pycache__

cleanall:
	rm -rf results
	rm -rf solutions
	rm -rf __pycache__
