import copy 
import numbers
import numpy as np
class Matrica():
    def __init__(self,rows=None,columns=None,epsilon=1e-9):
        self.epsilon = epsilon
        if rows and columns:
            self.matrix=[]
            for i in range(0,rows):
                temp = [0]*columns
                self.matrix.append(temp)
            self.init_rows_cols()
                    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if (np.allclose(np.array(self.matrix), np.array(other.matrix),atol=self.epsilon) and
                self.rows == other.rows and
                self.columns == other.columns):
                return True
        return False
    
    def __add__(self, other):
        if isinstance(other,Matrica):
            return self.add_matrix(other)
        elif isinstance(other, numbers.Number):
            return self.add_scalar(other)
        raise Exception("Objekt nije ni matrica ni skalar")
        
    def __radd__(self, other):
        return self.__add__(other)
    
    def __sub__(self, other):
        if isinstance(other,Matrica):
            return self.sub_matrix(other)
        elif isinstance(other, numbers.Number):
            return self.sub_scalar(other)
        raise Exception("Objekt nije ni matrica ni skalar")
        
    def __rsub__(self, other):
        return self.__sub__(other)
    
    def __mul__(self,other):
        if isinstance(other,Matrica):
            return self.multiply_matrices(other)
        elif isinstance(other, numbers.Number):
            return self.multiply_matrix_scalar(other)
        raise Exception("Objekt nije ni matrica ni skalar")
    
    def __rmul__(self, other): #right mul s matricama fix
        return self.__mul__(other)
    
    def __div__(self, other):
        if isinstance(other, numbers.Number):
            return self.divide_matrix_scalar(other)
        raise Exception("Objekt nije skalar")
    
    def __str__(self):
        longest_number=0
        for row in self.matrix:
            longest_number = max(longest_number,max(map(len,map(str,row)))) 
            #pretvorim sve elemente u string, pogledam njihovu duljinu i nadem duljinu max elementa 
            #i izaberem max od prijasnjeg poznatog i novog max elementa u retku
            
        string = ""
        for i in range(0,self.rows):
            string +="|"
            for j in range(0,self.columns):
                string +='%*s' % (longest_number+2,self.matrix[i][j])
            if i<self.rows-1:
                string +="|\n"
            else:
                string +="|"
        return string
                
    def load_matrix(self,file_path):
        with open(file_path,'r') as f:
            matrix=[]
            for line in f.readlines():
                row=[]
                row = line.strip().split()
                row = map(float,row)
                matrix.append(row)
        self.matrix = matrix
        self.init_rows_cols()
        
    def load_from_numpy(self,np_arr):
        matrix=[]
        for row in np_arr:
            matrix.append(row.copy())
        self.matrix = matrix
        self.init_rows_cols()
        
    
    def init_rows_cols(self):
        self.rows = len(self.matrix)
        self.columns = len(self.matrix[0])
        for row in self.matrix:
            if self.columns != len(row):
                raise Exception("Matrica nije pravilno zadana")
    
    def add_scalar(self,scalar):
        Mat = self.deep_copy()
        for i in range (0,self.rows):
            for j in range (0,self.columns):
                Mat.matrix[i][j] += scalar
        return Mat
    
    def sub_scalar(self,scalar):
        Mat = self.deep_copy()
        for i in range (0,self.rows):
            for j in range (0,self.columns):
                Mat.matrix[i][j] -= scalar
        return Mat
    
    def add_matrix(self,M):
        if M.rows != self.rows or M.columns != self.columns:
            raise Exception("Matrice nisu istih dimenzija")
        Mat = self.deep_copy()
        for i in range (0,self.rows):
            for j in range (0,self.columns):
                Mat.matrix[i][j] +=M.matrix[i][j]
        return Mat
        
    def sub_matrix(self,M):
        if M.rows != self.rows or M.columns != self.columns:
            raise Exception("Matrice nisu istih dimenzija")
        Mat = self.deep_copy()
        for i in range (0,self.rows):
            for j in range (0,self.columns):
                Mat.matrix[i][j] -=M.matrix[i][j]
        return Mat
    
    def multiply_matrix_scalar(self,scalar):
        Mat = self.deep_copy()
        for i in range (0,self.rows):
            for j in range (0,self.columns):
                Mat.matrix[i][j] *= scalar
        return Mat
    
    def divide_matrix_scalar(self, scalar):
        Mat = self.deep_copy()
        for i in range (0,self.rows):
            for j in range (0,self.columns):
                Mat.matrix[i][j] /= scalar
        return Mat
    
    def multiply_matrices(self,other):
        Mat = Matrica(self.rows,other.columns,epsilon=self.epsilon)
        if self.columns != other.rows:
            raise Exception("Matrices arent compatible for multiplication")
        for i in range(len(self.matrix)):
            for j in range(len(other.matrix[0])):
                for k in range(len(other.matrix)):
                    Mat.matrix[i][j] += self.matrix[i][k] * other.matrix[k][j]
        return Mat
    
    def deep_copy(self):
        matrix = []
        for i in range (0,self.rows):
            row = []
            for j in range (0,self.columns):
                row.append(self.matrix[i][j])
            matrix.append(row)
            
        M = Matrica(epsilon=self.epsilon)
        M.matrix = copy.deepcopy(matrix)
        M.init_rows_cols()
        return M
    
    def transpose(self):
        new_matrix = []
        for j in range(0,self.columns):
            row = []
            for i in range(0,self.rows):
                row.append(self.matrix[i][j])
            new_matrix.append(row)
        self.matrix = new_matrix
        self.init_rows_cols()
        return self
        
    
    
    def LU_decomposition(self):
        M = self.deep_copy()
        if self.rows != self.columns:
            raise Exception("Matrica nije kvadratna")
        for i in range (0,self.rows):
            for j in range (i+1, self.rows):
                if abs(M.matrix[i][i]) < self.epsilon:
                    raise Exception("Stozerni element je nula")
                M.matrix[j][i] /=1.* M.matrix[i][i]
                for k in range (i+1,self.rows):
                    M.matrix[j][k] -= M.matrix[j][i]*M.matrix[i][k]
        
        if abs(M.matrix[self.rows-1][self.rows-1]) < self.epsilon: #hacknuto treba sredit unutar petlje
            raise Exception("Stozerni element je nula")
        return M
                    
    def LUP_decomposition(self):
        M = self.deep_copy()
        P = []
        for i in range(0,self.rows):
            P.append(i)
        for i in range (0,self.rows):
            pivot = i;
            for j in range (i+1, self.rows):
                if abs(M.matrix[P[j]][i]) > abs(M.matrix[P[pivot]][i]):
                    pivot = j;
            if abs(M.matrix[P[pivot]][i]) < self.epsilon:
                raise Exception("Stozerni element je manji od zadane granice i dekompozicija je zaustavljena")
            temp = P[i]
            P[i] = P[pivot]
            P[pivot] = temp
            for j in range (i+1, self.rows):
                M.matrix[P[j]][i] /=1.* M.matrix[P[i]][i];
                for k in range (i+1, self.rows):
                    M.matrix[P[j]][k] -= M.matrix[P[j]][i] * M.matrix[P[i]][k];
        return M, P
    
    
    def return_L_sub_matrix(self):
        M = self.deep_copy()
        for i in range (0,self.rows):
            M.matrix[i][i]=1
            for j in range (i+1,self.rows):
                M.matrix[i][j]=0
        return M
    
    def return_U_sub_matrix(self):
        M = self.deep_copy()
        for i in range (0,M.rows-1):
            for j in range (i+1,M.rows):
                M.matrix[j][i]=0
        return M
    
    
    def forward_supstitution(self, system_matrix):
        if self.columns !=1:
            raise Exception("Solution matrix needs to have one column")
        if self.rows != system_matrix.rows:
            raise Exception("Solution vector doesn't have the same dimension as system matrix")
        M = self.deep_copy()
        for i in range (0,M.rows-1):
            for j in range (i+1, M.rows):
                M.matrix[j][0] -= system_matrix.matrix[j][i]*M.matrix[i][0]
        return M
        
    def backward_supstitution(self, system_matrix):
        if self.columns !=1:
            raise Exception("Solution matrix needs to have one column")
        if self.rows != system_matrix.rows:
            raise Exception("Solution vector doesn't have the same dimension as system matrix")
        M = self.deep_copy()
        for i in range (self.rows-1,-1,-1):
            M.matrix[i][0] /= system_matrix.matrix[i][i]
            for j in range (0, i):
                M.matrix[j][0] -= system_matrix.matrix[j][i]*M.matrix[i][0]
        return M
    
    def print_matrix(self):
        for i in range(0,self.rows):
            print(self.matrix[i])
        
    def arrange_in_order(self,P):
        if max(P)> self.rows-1:
            raise Exception ("Indexes of rows are too large")
        if min(P)< 0:
            raise Exception ("Indexes of rows are too small")
        if len(set(P)) != self.rows:
            raise Exception ("All indexes need to be inside the array P")
            
        matrix = []
        for index in P:
            matrix.append(copy.deepcopy(self.matrix[index]))
            
        M = self.deep_copy()
        M.matrix = matrix
        return M
    
    def return_as_numpy(self):
        return np.array(self.matrix).copy()
    
    def inverse(self):
        row = len(self.matrix)
        column = len(self.matrix[0])
        if row != column:
            raise Exception("Matrix doesnt have full rank")
        
        
        
        LUP_decomposed, P = self.LUP_decomposition()
        LUP_decomposed_arrange = LUP_decomposed.arrange_in_order(P)
        L = LUP_decomposed_arrange.return_L_sub_matrix()
        U = LUP_decomposed_arrange.return_U_sub_matrix()
        
        sol_array = np.array([])
        for i in range(row):
            sol = np.zeros((row,1))
            sol[i] = 1
            b = Matrica()
            b.load_from_numpy(sol)
            b = b.arrange_in_order(P)
            
            Y = b.forward_supstitution(L)
            X = Y.backward_supstitution(U)
            if not sol_array.any():
                sol_array = X.return_as_numpy().copy()
            else:
                sol_array = np.hstack((sol_array,X.return_as_numpy().copy()))
            
        return sol_array

def solve_matrix_by_LU_decomposition(system_matrix,solution_matrix):
    LU_decomposed = system_matrix.LU_decomposition()
    L = LU_decomposed.return_L_sub_matrix()
    U = LU_decomposed.return_U_sub_matrix()
    Y = solution_matrix.forward_supstitution(L)
    X = Y.backward_supstitution(U)
    return X

def solve_matrix_by_LUP_decomposition(system_matrix,solution_matrix):
    LUP_decomposed, P = system_matrix.LUP_decomposition()
    LUP_decomposed_arrange = LUP_decomposed.arrange_in_order(P)
    L = LUP_decomposed_arrange.return_L_sub_matrix()
    U = LUP_decomposed_arrange.return_U_sub_matrix()
    sol_matrix_arranged = solution_matrix.arrange_in_order(P)
    Y = sol_matrix_arranged.forward_supstitution(L)
    X = Y.backward_supstitution(U)
    return X
