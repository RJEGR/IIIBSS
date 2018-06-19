def addTip(bill):
    return 1.1 * bill

def equallyDivide(bill, n=3):
   total = addTip(bill)
   return total / n
