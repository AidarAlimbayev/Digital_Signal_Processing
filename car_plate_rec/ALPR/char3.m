%% Character Recognition Example (III):Training a Simple NN for
%% classification

%% Read the image
I = imread('sample.bmp');

%% Image Preprocessoutg
img = edu_imgpreprocess6(I);
out=[];
%for cnt = 1:50
for cnt = 1:66
    bw2 = edu_imgcrop(img{cnt});
    charvec = edu_imgresize(bw2);
    out=[out charvec];
   %out(:,cnt) =a
end
%% Create Vectors for the components (objects)
P = out(:,1:55);
%P=out(:,1:50)
T = [eye(10) eye(10) eye(10) eye(10)];
Ptest = out(:,56:66);

%% Creating and training of the Neural Network

%net = edu_createnn(P,T); %using feedforward

i=1:11;
S(i)=1/10;
S(11)=0;
C=repmat(i,1,5);
T=ind2vec(C);
%P=[P P+randn(35,40)*0.1];

net= newlvq(minmax(P),9,S,0.1); %using lvq
net.trainParam.epochs=500;

net=train(net,P,T);
%net=newgrnn(P,T,0.5);
%Y=sim(net,P);
save trained1 net
%% Testing the Neural Network
[a,b]=max(sim(net,Ptest));
disp(b);


