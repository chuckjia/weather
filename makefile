all:
	python weather.py

fresh: clean
	python weather.py

clean:
	rm -rf results
	rm -rf solutions
	rm -rf __pycache__

clear:
	clear
	clear
	clear
