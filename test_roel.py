base_array = ['a', 'b', 'c', 'd', 'e']
eval_array = np.array([0, 1, 1, 0, 1])

new = [x for (x, y) in zip(base_array, eval_array) if y == 1]

print new


new = [x for (x, y) in zip(base_array, eval_array) if y == 0]

print new

newBool = [x==0 for x in eval_array]

print newBool


for a, b, c in zip(base_array, eval_array, newBool):
    print a,b,c

for i in zip(base_array, eval_array, newBool):
    print 'dit is een tuple:', i


vb_tuple = ('a', 56845, ['dit', 'is', 'een', 'vector'])

a ,b, c = vb_tuple
print a, type(a)
print b, type(b)
print c, type(c)
