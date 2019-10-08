#Socring balance and its optimization
#Based on Kullback-Leibler
#Define x.log.y
x.log.y<-function(x,y){
	ifelse(x==0, 0, x*log(x/y, 2))
}
#KL for vectors x, y
divergence<-function(o, a){
	sum(x.log.y(a, o))
}
#Relative balance
r.balance<-function(o, a){
	1+divergence(o, a)/(log(min(o)[1], 2))
}
#Coverage
coverage<-function(No, Na, o, a){
	(Na/No)*r.balance(o, a)
}

#Optimization

#Table has columns:
#1. End groups
#2. Actual numbers
#3. Ideal percentages

#Argument n is the desired number of accessions to be allocated o removed


allocateAccessions<-function(table,n){
myset<-NA #Initialize set of end groups
e<-as.character(table[,1]) #end groups
x<-table[,2] #actual numbers
p<-table[,3]/100 #optimal proportions
	for(i in 1:n){
	a<-NA #initialize vector of relative balance scores
		for(j in 1:length(x))
		{b<-x;b[j]<-x[j]+1;a[j]<-r.balance(p,b/sum(b))}
		index<-which(a==max(a))[1]
		myset[i]<-e[index]
		x[index]<-x[index]+1
	}
a1<-table(myset)
a1<-as.data.frame(a1)
names(a1)<-c("end.group","allocate")
a2<-r.balance(p,x/sum(x))
list("allocation" = a1,"RB" = a2)
}

removeAccessions<-function(table,n){
myset<-NA #Initialize set of end groups
e<-as.character(table[,1]) #end groups
x<-table[,2] #actual numbers
p<-table[,3]/100 #optimal proportions
	for(i in 1:n){
	a<-NA #initialize vector of relative balance scores
		for(j in 1:length(x))
		{b<-x;b[j]<-x[j]-1;a[j]<-r.balance(p,b/sum(b))}
		index<-which(a==max(a))[1]
		myset[i]<-e[index]
		x[index]<-x[index]-1
	}
a1<-table(myset)
a1<-as.data.frame(a1)
names(a1)<-c("end.group","remove")
a2<-r.balance(p,x/sum(x))
list("allocation" = a1,"RB" = a2)
}


#This is not optimization, but simple math. 
#Decide a target total number of accessions and do the actions to balance the whole collection.
#The resulting collection will be perfectly balanced
#n is the target collection size

makeEven<-function(table,n){
ideal.n<-n*table[,3]/100
dif<-ideal.n-table[,2]
dif<-round(dif,0)
data.frame("end"=table[,1],"action"=dif)
}




