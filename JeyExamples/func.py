class structure():

    def __init__(self):
        self.x=46
        self.y=23
        self.z=25
        print("hello _init_")

    def initalise(self):
        self.c = self.x * self.y + self.z
    
    def output(self):
        print(self.c)

if __name__ == "__main__":
    print("main call !")
    sim = structure()

print("jeyaramanan")