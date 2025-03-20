# importing packages 
import scipy.stats as stats

# creating data 
data = [[252448, 1500], [162071402.3, 32184.00898]]

# performing fishers exact test on the data 
odd_ratio, p_value = stats.fisher_exact(data)
print('odd ratio is : ' + str(odd_ratio))
print('p_value is : ' + str(p_value)) 