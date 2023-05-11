import multiprocessing as mp
mp.set_start_method('fork')

def sayhello():
    print("Hello!")
    
    
job = mp.Process(target=sayhello)
job.start()

# Wait for job to finish
job.join()

print("All done ")
print(mp.get_start_method())