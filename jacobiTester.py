import os
import time
import pandas as pd

parse = lambda x: "".join(x.split()) # Remove Whitespaces.
PYTHON_IMPL = "python3 spkmeans.py"
C_IMPL = "./spkmeans"
txt = ".txt"
csv = ".csv"
input = " 0 jacobi Jacobi_tests/matrix"
output = "Jacobi_tests/jacobi_output"
list = list(map(str, range(2,51)))
tests = [input + index + csv for index in list]
outputs = [output+ index + txt for index in list]
count = 2
df = pd.DataFrame() 

# c
for output, test in zip(outputs, tests):

	c_start = time.time()
	c_output = os.popen(f"{C_IMPL}{test}").read()
	c_end = time.time()
	bool = ' '
	with open(output, "r") as f:
		my_output = f.read()

	if parse(c_output) != parse(my_output):
		bool = 'F'

	dict = {'C run time':c_end - c_start,'C result':bool,'inputs':f'matrix{count}'}
	df = df.append(dict,ignore_index=True)
	count = count+1

count = 2
p_runtime = []
p_result = []
# python
for output, test in zip(outputs, tests):

	python_start = time.time()
	python_output = os.popen(f"{PYTHON_IMPL}{test}").read()
	python_end = time.time()
	bool = ' '

	with open(output, "r") as f:
		my_output = f.read()

	if parse(python_output) != parse(my_output):
		bool = 'F'

	p_runtime.append(python_end - python_start)
	p_result.append(bool)
	count = count+1

df['python result'] = p_result
df['python run time'] = p_runtime
df = df.reindex(columns=['inputs','C result','python result','C run time','python run time'])
pd.DataFrame(df).to_csv('jacobi tester result.csv', index=False)
print('End')

