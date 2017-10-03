class Tester(dict):
    def __init__(self, items):
        super(Tester, self).__init__()


a = dict(zip(('a', 'b', 'c'), range(1, 4)))

print(a)

b = Tester(a)

print(b)