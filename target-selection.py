#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:25:17 2020

@author: lizeth.machado
"""
import os
import glob
import nibabel as nib
import numpy as np
from skimage import measure
import math
from scipy import signal
import argparse
#Copy the python program to your files from: 
parser = argparse.ArgumentParser()
parser.add_argument( '--target_dir',
                    nargs='?',
                    required=True,
                    help='path to the dir that contains images, it will look for all .nii an d.nii.gz images')

args = parser.parse_args()
x=args.target_dir
os.chdir(x)

#-----------Variables 

images=[] 
first=[]
last=[]
middle=[]
Good=[]
Good_lst=[]
list_blackSlides=[]
list_mot_ssim=[]
list_mot_mse=[]

    #--- Parameter 2------
sigma=[]
Gauss_avg=[]
dice_avg=[]
correlation_avg=[]

M = [[1, -2, 1],
   [-2, 4, -2],
   [1, -2, 1]]

#get name in the current directory
get_files=glob.glob('*.nii')
get_files_gz=glob.glob('*.nii.gz')
get_files += get_files_gz


ima_1=[ x for x in get_files if "mm.nii" in x ]
ima_1.sort()
for n in ima_1:
    #--- Parameter 1------
    flag_ac1=False
    flag_ac2=False
    block=False
    count=[]
    mse=[]
    ssim=[]
    o=0
    j=0
    #--- Parameter 2------
    dice=[]
    cor=[]
    image2d=[]
    result=[]
    sigma=[]
    count1=0
    l=0
    c=0
    
    size=[nib.load(n).get_data()][0].shape[2]
    for i in range(size):
        img1=np.asarray([nib.load(n).get_data()[:,:,i]])[0,:,:]
        images.extend(img1)
        if(np.count_nonzero(img1)!=0):
            flag_ac=False
            count.append(np.count_nonzero(img1))
            o+=1
        else:
            flag_ac=True
            
        if flag_ac1==True and flag_ac== False:
            flag_ac2=True 
            
        if flag_ac2 ==True and block==False:
            j=i;
            block=True;
            
        flag_ac1=flag_ac
        
        #-----------------------Parameters to order array-----------
        if i<(size-1):
            img2=np.asarray([nib.load(n).get_data()[:,:,i+1]])[0,:,:]
            ssim.append(measure.compare_ssim(img1, img2))
            mse.append(measure.compare_nrmse(img1, img2))
        
        H,W=np.shape(img1)
        nonzero=np.count_nonzero(img1)
        
        if nonzero==0:
            l+=1
            
        if nonzero!=0:
            temp=(np.sum(np.sum(np.absolute(signal.convolve2d(img1, M)))))
            sigma.append(temp * math.sqrt(0.5 * math.pi) / (6 * (W-2) * (H-2)))

            if count1<7:
                h=int((size-l)/2-3)
                img1=np.asarray([nib.load(n).get_data()[:,:,(i+h)]])
                img1=img1[0,:,:]
                im1 = np.asarray(img1).astype(np.bool)
                img12=np.asarray([nib.load(n).get_data()[:,:,(i+h)+1]])
                img12=img12[0,:,:]
                im2 = np.asarray(img12).astype(np.bool)
                im_sum = im1.sum() + im2.sum()
                # Compute Dice coefficient
                intersection = np.logical_and(im1, im2)
                dice.append(2. * intersection.sum() / im_sum)
                #----------------Correlation Ratio------------------
                x=img1
                y=img12
                mu_x = x.mean(1)
                mu_y = y.mean(1)
                dof = x.shape[1]
                s_x = x.std(1, ddof=dof - 1)
                s_y = y.std(1, ddof=dof - 1)
                cov = np.dot(x,y.T) - dof * np.dot(mu_x[:, np.newaxis],mu_y[np.newaxis, :])
                sqrt=np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])
                result.append(np.sum(cov)/np.sum(sqrt))
                count1+=1

        
    correlation_avg.extend(sum(result)/np.shape(result))
    dice_avg.append(sum(dice)/7.0)         
    Gauss_avg.extend(sum(sigma)/np.shape(sigma))   
    #Check if there are black slides in the middle of the brain 
    y=np.argwhere(np.isnan(mse))[:,0]
    y=[x - y for x, y in zip(y[1:],y)]
    y=np.count_nonzero([x>1 and x<45 for x in y]) #45 because it can have black slides in the beggininand at the end
    #check motion 
    mot_ssim=np.count_nonzero([x<0.91 for x in ssim])
    mot_mse=np.count_nonzero([x>2 and x<1.8E308 for x in mse])
     #list for all the images 
    list_blackSlides.append(y)
    list_mot_ssim.append(mot_ssim)
    list_mot_mse.append(mot_mse)
            
    first.append(sum(count[:2])/2)
    last.append(sum(count[len(count)-2:])/2)
    ima=np.asarray([nib.load(n).get_data()[:,:,(int(j+(o/2)))]])
    ima=ima[0,:,:]
    middle.append(np.count_nonzero(ima))


#------------------Results----------------------    
    #Complete brains
ratio_fst=np.divide(first,middle)
ratio_lst=np.divide(last,middle)
Good=[index for index,value in enumerate(ratio_fst) if value <0.35]
Good_lst=[index for index,value in enumerate(ratio_lst) if value <0.35]

Good=np.intersect1d(Good,Good_lst)
    #Get the order of the images 
index_bs=np.argsort(list_blackSlides)# excludes
index_mot_ssim=np.argsort(list_mot_ssim)
index_mot_mse=np.argsort(list_mot_mse) 
    #Parameters for best images
arr=np.array(correlation_avg)
idx3=(-arr).argsort()[:3]

arr=np.array(Gauss_avg)
idx2=(-arr).argsort()[np.shape(Gauss_avg)[0]-3:]

arr=np.array(dice_avg)
idx1=(-arr).argsort()[:3]

sim1=np.intersect1d(idx1,idx2)
sim2=np.intersect1d(idx2,idx3)
sim3=np.intersect1d(idx1,idx3)
inter=np.intersect1d(index_mot_mse[:np.int(len(index_mot_mse)/2)],Good)
inter=np.intersect1d(index_mot_ssim[:np.int(len(index_mot_ssim)/2)],inter)
inter = [i for i in inter if i not in index_bs[np.int(len(index_bs))-3:]]
result=np.intersect1d(inter,sim3)


print('----------------------------Result------------------')
os.system("sed -n '14p' run_reconstruction.sh")
if len(result)!=0:
    print(result)
    print(ima_1[result[0]])
else:
    result=np.intersect1d(inter,sim1)
    print(inter,sim3,sim1)
    if len(result)!=0:
        print('The result is: '+str(result[0]))
        print(ima_1[result[0]])
    else: 
        result=np.intersect1d(inter,sim2)
        print(inter,sim3,sim1,sim2)
        if len(result)!=0:
            print('The result is: '+str(result[0]))
            print(ima_1[result[0]])
        else:
            print('THERE IS NO MATCH!')
    
