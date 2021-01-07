
run:
	python3 app.py compute --name sample --error-dist 2 &

server-run:
	nohup python3 app.py compute --name sample --error-dist 2 &

test:
	pytest