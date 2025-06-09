CSN<-function(data){
  weighted = 0
  boxsize = 0.1
  alpha = 0.01
  n1=nrow(data)
  n2=ncol(data)
  c=1:n2
  upper<-matrix(0,n1,n2,dimnames = list(row.names(data),
                                        colnames(data)))
  lower<-matrix(0,n1,n2,dimnames = list(row.names(data),
                                          colnames(data)))
  for (i in 1:n1) {
    s1<-sort(data[i,])
    s2<-order(data[i,])
    s3<-ifelse(s1>0,1,ifelse(s1<0,-1,0))
    n3=n2-sum(s3)
    h=round(boxsize/2*sum(s3))
    k=1
    while(k<=n2){
      s=0
      while(k+s+1 <= n2 && s1[k+s+1] == s1[k]){s=s+1}
      if(s>=h){
        upper[i,s2[k:(k+s)]]=data[i,s2[k]]
        lower[i,s2[k:(k+s)]]=data[i,s2[k]]
      }else{
        upper[i,s2[k:(k+s)]]=data[i,s2[min(n2,(k+s+h))]]
        lower[i,s2[k:(k+s)]]=data[i,s2[max(n3*(n3>h)+1,(k-h))]]
      }
      k=k+s+1
    }
  }
  csn<-list()
  p=2.3263
  B<-matrix(0,n1,n2,dimnames = list(row.names(data),
                                    colnames(data)))
  for(k in c){
    for (j in 1:n2) {
      B[,j]=data[,j]<=upper[,k]&data[,j]>=lower[,k]
    }
    a=apply(B,1,sum)
    eps=2.2204e-16
    d=(B%*%t(B)*n2-a%*%t(a))/sqrt((a%*%t(a))*((n2-a)%*%t((n2-a)))/(n2-1)+eps)
    for (i in 1:nrow(d)) {
      d[i,i]=0
    }
    csn[[k]]<-d*(d>p)
    print(paste("Cell",k,"is completed"))
  }
  csn
}
