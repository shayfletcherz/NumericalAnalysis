#----------------------
#Shay Fletcher     318727641
#Nika Tatikishvili 321735433
#----------------------

#פונקציה אשר יוצרת מטריצה ריבועית לפי גודל מסוים של מטריצת יחידה
def identityMat(n):
    mat = [([0] * n) for i in range(n)]  #יצירת מטריצת אפסים
    for i in range(0, n):
        mat[i][i] = 1  #הפיכת אפסים ל-1 באלכסון
    return mat

#בדיקה האם פונקציה סינגולרית בעזרת חישוב דטרמיננטה
def detCalc(matrix):
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    det = 0
    for j in range(0, n):
        det += ((-1) ** j) * matrix[0][j] * detCalc(detReduce(matrix, 0, j))
    return det

#הפונקציה תיצור מטריצה שמחסרת את העמודות והשורות שלה על מנת לעזור בחישוב הדטרמיננטה
def detReduce(matrix, i, j):
    return [row[:j] + row[j + 1:] for row in (matrix[:i] + matrix[i + 1:])]


def multiplyMatrix(mat1, mat2):
    if len(mat1[0]) != len(mat2):
        print("Illegal multiplication")
    finalMat = [([0] * len(mat2[0])) for i in range(len(mat1))]

    #מעבר על שורות מטריצה במטריצה הראשונה
    for row1 in range(0, len(mat1)):
        # מעבר על עמודות במטריצה השנייה
        for col2 in range(0, len(mat2[0])):
            # מעבר על שורות המטריצה השנייה
            for row2 in range(0, len(mat2)):
                finalMat[row1][col2] += mat1[row1][row2] * mat2[row2][col2]
    return finalMat

#חישוב חיבור מטריצות
def addMatrix(mat1, mat2):
    size = len(mat1)
    finalMat = [([0] * size) for i in range(size)]

    for row in range(0, len(mat1)):
        for col in range(0, len(mat1[0])):
            finalMat[row][col] = mat1[row][col] + mat2[row][col]
    return finalMat

#יצירת העתק למטריצה המקורית לצורך חישובים
def copyMat(mat):
    size = len(mat)
    copyMat = [([0] * size) for i in range(size)]
    for row in range(0, len(mat)):
        for col in range(0, len(mat[0])):
            copyMat[row][col] = mat[row][col]
    return copyMat

#פונקציה שיוצרת מטריצת אלמנטרית ובוחר מיקום מסוים על מנת לאפס אותו
def reset(row, col, x, pivot, zeroval):
    #תחילה המטריצה תוצג כמטריצת יחידה
    #x מייצג את סך הכל השורות ועמודות
    elementaryMat = identityMat(x)
    #במיקום הנבחר באותה מטריצה נמקם ערך שיציג לנו 0 במידה ונכפיל בין המטריצות
    #zeroval תייצג את המקום שנרצה לאפס
    #pivot מייצג את הפיבוט עליו המטריצה משתמשת
    elementaryMat[row][col] = -(zeroval / pivot)
    return elementaryMat

#פונקציית עזר שתאפשר לנו ליצור מטריצה אלמנטרית שתאפשר לנו להכפיל שורה לפי ערך הפיבוט
def multiplyRow(row, val, x):
    elementaryMatrix = identityMat(x)
    elementaryMatrix[row][row] = val
    return elementaryMatrix

#הפונקציה תייצר מטריצה אלמנטרית שבעזרתה נוכל לחליף בין השורות בעזרת כפל בהם
def swapRow(row1, row2, x):

    elementaryMatrix = identityMat(x)
    elementaryMatrix[row1][row1] = 0
    elementaryMatrix[row1][row2] = 1
    elementaryMatrix[row2][row2] = 0
    elementaryMatrix[row2][row1] = 1
    return elementaryMatrix

#פונקציה להדפסת מטריצה הופכית
def InvertMatrix(matrix):
    if len(matrix) != len(matrix[0]):
        print("singular matrix!")
    n = len(matrix)
    inverted = identityMat(n)
    for j in range(0, n):
        for i in range(0, n):
            if i == j:
                pivot = matrix[i][j]
                for k in range(i + 1, n):
                    if abs(matrix[k][j]) > abs(pivot):  #ביצוע הפיבוט
                        elementaryMatrix = swapRow(k, i, n)
                        matrix = multiplyMatrix(elementaryMatrix, matrix)
                        inverted = multiplyMatrix(elementaryMatrix, inverted)
                        pivot = matrix[i][j]

                if matrix[i][j] == 0:
                    print("singular matrix!")

        for i in range(0, n):
            if i != j:
                if matrix[i][j] != 0:
                    elementaryMatrix = reset(i, j, n, pivot, matrix[i][j])
                    matrix = multiplyMatrix(elementaryMatrix, matrix)
                    inverted = multiplyMatrix(elementaryMatrix, inverted)

    for i in range(0, n):
        if matrix[i][i] != 1:
            if matrix[i][i] < 0:
                elementaryMatrix = multiplyRow(i, -1, n)
                matrix = multiplyMatrix(elementaryMatrix, matrix)
                inverted = multiplyMatrix(elementaryMatrix, inverted)

            elementaryMatrix = multiplyRow(i, 1 / matrix[i][i], n)
            matrix = multiplyMatrix(elementaryMatrix, matrix)
            inverted = multiplyMatrix(elementaryMatrix, inverted)
    for row in range(n):
        for col in range(n):
            inverted[row][col] = round(inverted[row][col], 2)
    return inverted

#חישוב LU בעזרת L וU
def luMatrixCalculation(mat):

    size = len(mat)
    res = copyMat(mat)
    for col in range(0, size):
        pivot = res[col][col]
        if pivot == 0 and col != size - 1:
            for row in range(col + 1, size):
                if res[row][col] != 0:
                    elementary_matrix = swapRow(row, col, size)
                    res = multiplyMatrix(elementary_matrix, res)
                    mat = multiplyMatrix(elementary_matrix, mat)
                    pivot = res[col][col]
                    break
        if pivot == 0 and col != size - 1:
            print("Can not calculate!")

        for row in range(col + 1, size):
            m = -(res[row][col] / pivot) #מספר כופל
            elementary_mat = identityMat(size) #יצירת מטריצה יחידה בגודל המטריצה המקורית
            elementary_mat[row][col] = m
            res = multiplyMatrix(elementary_mat, res) #בסוף הריצה בלולאה נקבל את מטריצה U
            if col == 0 and row == 1:
                L = InvertMatrix(elementary_mat)
            else:
                L = multiplyMatrix(L, InvertMatrix(elementary_mat))

    if multiplyMatrix(L, res) == mat:
        print("\nL Matrix:", L)
        print("U Matrix:", res)
    else:
        raise Exception("Can not calculate!")

#מציאת וקטור תוצאה
def findResultVector(matrix, b):
    print("\nThe result vector is: ")
    return multiplyMatrix(InvertMatrix(matrix), b)


#בדיקה האם מטריצה סינגולרית או לא סינגולרית ומחשבת LU או מוצאת וקטור תוצאה
def CalcMatrix(matrix, b):
    if detCalc(matrix) == 0:
        luMatrixCalculation(matrix)
    else:
        print(findResultVector(matrix, b))
        #print("Cond(A)=", end=' ')
        #print(CalcConditionNumber(matrix))

#מטריצה הפיכה - לפי גאוס
mat1 = [[1, -1, -2],
        [2, -3, -5],
        [-1, 3, 5]]

b = [[1], [2], [3]]

#מטריצה לא הפיכה - LU
mat2 = [[1, 1, 5],
        [1, 2, 7],
        [2, -1, 4]]

#driver
try:
    CalcMatrix(mat1, b)
    CalcMatrix(mat2, b)
except Exception as e:
    print(e)
