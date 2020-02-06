#!/usr/bin/env python                                            
#
# mareconal ds ChRIS plugin app
#
# (c) 2016-2019 Fetal-Neonatal Neuroimaging & Developmental Science Center
#                   Boston Children's Hospital
#
#              http://childrenshospital.org/FNNDSC/
#                        dev@babyMRI.org
#



"-----------------Libaries---------------"
import os,sys
import glob
import nibabel as nib
import matplotlib.pyplot as plt
import argparse
#import numpy as np
from PIL import Image 
import math
import numpy as np
from numpy import zeros
import pudb

    
sys.path.append(os.path.dirname(__file__))



class Mareconal():
    """
    This app does masking, reconstruction and alignment  .
    """
    AUTHORS                 = 'Lizeth (lizmachado09@gmail.com)'
    SELFPATH                = os.path.dirname(os.path.abspath(__file__))
    SELFEXEC                = os.path.basename(__file__)
    EXECSHELL               = 'python3'
    TITLE                   = 'Mareconal: Part of the fetal reconstruction pipeline'
    CATEGORY                = ''
    TYPE                    = 'ds'
    DESCRIPTION             = 'This app does masking, reconstruction and alignment  '
    DOCUMENTATION           = 'http://wiki'
    VERSION                 = '0.1'
    ICON                    = '' # url of an icon image
    LICENSE                 = 'Opensource (MIT)'
    MAX_NUMBER_OF_WORKERS   = 1  # Override with integer value
    MIN_NUMBER_OF_WORKERS   = 1  # Override with integer value
    MAX_CPU_LIMIT           = '' # Override with millicore value as string, e.g. '2000m'
    MIN_CPU_LIMIT           = '' # Override with millicore value as string, e.g. '2000m'
    MAX_MEMORY_LIMIT        = '' # Override with string, e.g. '1Gi', '2000Mi'
    MIN_MEMORY_LIMIT        = '' # Override with string, e.g. '1Gi', '2000Mi'
    MIN_GPU_LIMIT           = 0  # Override with the minimum number of GPUs, as an integer, for your plugin
    MAX_GPU_LIMIT           = 0  # Override with the maximum number of GPUs, as an integer, for your plugin

    # Use this dictionary structure to provide key-value output descriptive information
    # that may be useful for the next downstream plugin. For example:
    #
    # {
    #   "finalOutputFile":  "final/file.out",
    #   "viewer":           "genericTextViewer",
    # }
    #
    # The above dictionary is saved when plugin is called with a ``--saveoutputmeta``
    # flag. Note also that all file paths are relative to the system specified
    # output directory.
    OUTPUT_META_DICT = {}

    "--------------Parse Arguments------------------"
    parser = argparse.ArgumentParser()
    parser.add_argument( '--target_dir',
                        nargs='?',
                        required=True,
                        help='path to the dir that contains images, it will recursivelly look for all .nii images and starts to run the pipeline')
    parser.add_argument( "--one", 
                        action="store_true",default= False ,
                        help= "if the option is specified you will only run one function and not the whole pipeline ")
    parser.add_argument("-f","--function", 
                        action="store",
                        help="Type the function you wan to run as 'function()' \r \
                        \nThe avaiable options are: \r \
                        \nverification_image() \r\
                        \nverification_image_mask() \r\
                        \nverification_image_mask_traspond() \r\
                        \nbias_field_correction() \r\
                        \nwriting_script()\r\
                        \nalignment()")
    args = parser.parse_args()
    target_dir = args.target_dir
    function=args.function
    os.chdir(target_dir)
    
    "--------------Functions------------------"
    
    #     
    def get_images():
        get_files=glob.glob('*.nii')
        get_files_gz=glob.glob('*.nii.gz')
        get_files += get_files_gz
                
        ima=[ x for x in get_files if "brain" not in x ]
        ima=[ x for x in ima if "mask" not in x ]
        ima=[ x for x in ima if "mm.nii" not in x ]
        ima=[ x for x in ima if "recon" not in x ]
        ima=[ x for x in ima if "corrected" not in x ]
        ima.sort()
        return ima
    def veri_nii_gz():
        flag=0
        get_files_gz=glob.glob('*.nii.gz')
        if len(get_files_gz)!=0:
            flag=1
        return flag
        
    def get_images_mask():
        get_files=glob.glob('*.nii')
        get_files_gz=glob.glob('*.nii.gz')
        get_files += get_files_gz
                
        ima=[ x for x in get_files if "mask" in x ]
        ima.sort()
        return ima
    
    def verification_image():
        # Variables 
        ima=get_images()
        images=[]
        #Load images
        for i in ima:
             size=[nib.load(i).get_data()][0].shape[2]
             images.extend([nib.load(i).get_data()[:,:,int(size/3)]])
             images.extend([nib.load(i).get_data()[:,:,int(size/2)]])
             images.extend([nib.load(i).get_data()[:,:,int(size-5)]])
        print('Creating images')   
         #Set properties of images  
        plt.rcParams['figure.facecolor'] = 'black'
        fig=plt.figure(figsize=(25, 15))
        rows =int(round((len(images)/3)/2.0))
        columns = int(round(len(images)/rows)+1)
        for i in range(1,len(images)+1):
            fig.add_subplot(rows, columns, i)
            plt.imshow(images[i-1], cmap='gray',interpolation='nearest')
            plt.axis('off')
        plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
    
        fig.show()
        plt.savefig('Im.png', dpi=300)
    
    def brain_masking_tool():
        #--------------Command from shell-------------
        #Getting directories
        cwd = os.getcwd()
        home=os.getenv("HOME")
        os.chdir(home)
        os.system('./brain_mask --target-dir ' +cwd)
        os.chdir(cwd)
        
    #----------------------MASK-----------------------------
    def verification_image_mask():
        images=[]
        name_images_mask=get_images_mask()
        #Load images
        for i in name_images_mask:
            size=[nib.load(i).get_data()][0].shape[2]
            images.extend([nib.load(i).get_data()[:,:,int(size/3)]])
            images.extend([nib.load(i).get_data()[:,:,int(size/2)]])
            images.extend([nib.load(i).get_data()[:,:,size-5]])
        print('Creating images')
        #Set properties of images
        plt.rcParams['figure.facecolor'] = 'black'
        fig=plt.figure(figsize=(25, 15))
        rows =int(round((len(images)/3)/2.0))
        columns = int(round(len(images)/rows)+1)
        for i in range(1,len(images)+1):
            fig.add_subplot(rows, columns, i)
            plt.imshow(images[i-1], cmap='gray',interpolation='nearest')
            plt.axis('off')
          
        plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
        fig.show()
        plt.savefig('Im_mask.png', dpi=300)
    #------------------------resize---------------------- 
        # Opens a image 
        im = Image.open("Im_mask.png")  
          
        # Size of the image in pixels (size of orginal image)  
        width, height = im.size  
              
        newsize = (7500, 4500) 
        im1 = im.resize(newsize) 
        # Shows the image in image viewer  
        #im1.show()
        im1 = im1.save("Im_mask1.png")
    
    def verification_image_mask_traspond():
    #--------------------Traspond----------------------
        overlay = Image.open('Im.png')
        base = Image.open('Im_mask.png')
        
        bands = list(overlay.split())
        if len(bands) == 4:
            # Assuming alpha is the last band
            bands[3] = bands[3].point(lambda x: x*0.6)
        overlay = Image.merge(overlay.mode, bands)
        
        base.paste(overlay, (0, 0), overlay)
        base.save('result.png')
        
    def checking_mask():
        ima=get_images()
        os.system('. neuro-fs stable')
        for x in  range(len(ima)):
            os.system('~/arch/Linux64/packages/itksnap-2.4.0-20121121-Linux-x86_64/bin/itksnap '+ima[x])
            
    def bias_field_correction():
        ima=get_images()
        ima_mask=get_images_mask()
        flag=veri_nii_gz()
        for x in  range(len(ima)):
            try:
                if flag==1:
                    split=len(ima[x])-3
                else:
                    split=len(ima[x])
                os_cmd='mri_mask '+ima[x]+" "+ima_mask[x]+' brain_'+ima[x][:split]
                if os.system(os_cmd) != 0:
                    raise Exception('mri_mask not working')
                os.system('~/arch/Linux64/packages/ANTs/current/bin/N4BiasFieldCorrection -d 3 -o corrected_'+ima[x][:split]+' -i brain_'+ima[x][:split]+' -s 3  -c [400x400x400,0.00]')
            except:
                print("Oops!",sys.exc_info()[0],"occured.")
                print("AN ERROR HAS OCCURED: mri_mask tool not runing\n Freesurfer environment is not loaded")
                
    
    
    def writing_script():
        #-----------------Exctrating information--------------
        get_files=glob.glob('*.nii') 
        ima=[ x for x in get_files if "corrected" in x ]
        cwd = os.getcwd()
        string=""
        thickness=' '
        packages=''
        idstr=''
        imageDir=""
        imageDir1=""
        for x in  range(len(ima)):
            try:
                output=os.popen("mri_info "+ima[x]+" | sed -n 's/.*voxel sizes: //p'").read()
                thick=str(output[len(output)-7:len(output)-4])
                thickness+=str(output[len(output)-7:len(output)-4])
                split=len(ima[x])-4
                os.rename(ima[x],ima[x][10:split]+'-'+thick+'mm'+ima[x][split:]) 
                thickness+=' '
                packages+='1 '
                idstr+='id '
                imageDir+=" $imageDir/"+ima[x][10:split]+'-'+thick+'mm'+ima[x][split:]
                imageDir1+=" "+cwd+"/"+ima[x][10:split]+'-'+thick+'mm'+ima[x][split:]
            except:
                print("Oops!",sys.exc_info()[0],"occured.")
        string="/neuro/arch/Linux64/packages/irtk/build-tbb-rc/bin/reconstruction recon.nii "+str(len(ima))+imageDir1+" "+idstr+"-thickness "+thickness+" -debug -packages "+packages
        print(string)
    #--------------Writing script ---------------
        f= open("run_reconstruction.sh","w+")
        f.write("#!/bin/bash \r\
                \nappDirTBB='/neuro/arch/Linux64/packages/irtk/build-tbb-rc/bin\' \
                \nimageDir="+cwd+" \r \
                \ntransformDir='' \r \
                \nmkdir recon \
                \ncd recon \r \
                \n$appDirTBB/reconstruction recon.nii "+str(len(ima))+" "+imageDir+
                " "+idstr+"-thickness "+thickness+" -debug -packages "+packages)
        f.close()
        
        #writing the commands in terminal
        os.system("transformDir='' ")
        if not os.path.exists('./recon'):
            os.makedirs('recon')
    
        os.chdir('recon')
        os.system(string)
        
    def alignment():
        os.chdir('..')     
        path=os.getcwd()
        if not os.path.exists('./alignment'):
            os.makedirs('alignment')
            
        if not os.path.exists('./segmentation'):
            os.makedirs('segmentation')
            
        os.system('cp recon/recon.nii alignment')
        os.chdir(os.getenv("HOME"))
        
        if not os.path.exists('./templates'):
            os.system('cp -R /neuro/labs/grantlab/users/maria.hernandez/Trial/templates ~')
        
        os.chdir(path)
        os.chdir('alignment')
        
        os.system(\
        'flirt -in ~/templates/template-23/template-23.nii -ref recon.nii \
        -out Temp-Recon-7dof-23.nii -omat Temp-Recon-7dof-23.xfm -dof 7; \
        flirt -in ~/templates/template-23/csf-23.nii -ref recon.nii \
        -out csf-aligned23.nii -init Temp-Recon-7dof-23.xfm -applyxfm; \
        flirt -in ~/templates/template-24/template-24.nii -ref recon.nii \
        -out Temp-Recon-7dof-24.nii -omat Temp-Recon-7dof-24.xfm -dof 7; \
        flirt -in ~/templates/template-24/csf-24.nii -ref recon.nii \
        -out csf-aligned24.nii -init Temp-Recon-7dof-24.xfm -applyxfm; \
        flirt -in ~/templates/template-25/template-25.nii -ref recon.nii \
        -out Temp-Recon-7dof-25.nii -omat Temp-Recon-7dof-25.xfm -dof 7; \
        flirt -in ~/templates/template-25/csf-25.nii -ref recon.nii \
        -out csf-aligned25.nii -init Temp-Recon-7dof-25.xfm -applyxfm; \
        flirt -in ~/templates/template-26/template-26.nii -ref recon.nii \
        -out Temp-Recon-7dof-26.nii -omat Temp-Recon-7dof-26.xfm -dof 7; \
        flirt -in ~/templates/template-26/csf-26.nii -ref recon.nii \
        -out csf-aligned26.nii -init Temp-Recon-7dof-26.xfm -applyxfm; \
        flirt -in ~/templates/template-27/template-27.nii -ref recon.nii \
        -out Temp-Recon-7dof-27.nii -omat Temp-Recon-7dof-27.xfm -dof 7; \
        flirt -in ~/templates/template-27/csf-27.nii -ref recon.nii \
        -out csf-aligned27.nii -init Temp-Recon-7dof-27.xfm -applyxfm; \
        flirt -in ~/templates/template-28/template-28.nii -ref recon.nii \
        -out Temp-Recon-7dof-28.nii -omat Temp-Recon-7dof-28.xfm -dof 7; \
        flirt -in ~/templates/template-28/csf-28.nii -ref recon.nii \
        -out csf-aligned28.nii -init Temp-Recon-7dof-28.xfm -applyxfm; \
        flirt -in ~/templates/template-29/template-29.nii -ref recon.nii \
        -out Temp-Recon-7dof-29.nii -omat Temp-Recon-7dof-29.xfm -dof 7; \
        flirt -in ~/templates/template-29/csf-29.nii -ref recon.nii \
        -out csf-aligned29.nii -init Temp-Recon-7dof-29.xfm -applyxfm; \
        flirt -in ~/templates/template-30/template-30.nii -ref recon.nii \
        -out Temp-Recon-7dof-30.nii -omat Temp-Recon-7dof-30.xfm -dof 7; \
        flirt -in ~/templates/template-30/csf-30.nii -ref recon.nii \
        -out csf-aligned30.nii -init Temp-Recon-7dof-30.xfm -applyxfm; \
        flirt -in ~/templates/template-31/template-31.nii -ref recon.nii \
        -out Temp-Recon-7dof-31.nii -omat Temp-Recon-7dof-31.xfm -dof 7; \
        flirt -in ~/templates/template-31/csf-31.nii -ref recon.nii \
        -out csf-aligned31.nii -init Temp-Recon-7dof-31.xfm -applyxfm; \
        flirt -in ~/templates/template-32/template-32.nii -ref recon.nii \
        -out Temp-Recon-7dof-32.nii -omat Temp-Recon-7dof-32.xfm -dof 7; \
        flirt -in ~/templates/template-32/csf-32.nii -ref recon.nii \
        -out csf-aligned32.nii -init Temp-Recon-7dof-32.xfm -applyxfm;')
        
        recon = nib.load('recon.nii') # Load reconstruction image 
        size=recon.get_fdata().shape # Get the dimensions of the volume 
        
        meas = [0,0,0]
        beginning=[0,0,0]
        
        for i in range (0, 3):
            beginning[i]=int(round(size[i]/10.0))
            temp=size[i]-beginning[i]
            meas[i]=temp-beginning[i]
            
            if (i==0):
                dimensions=zeros([size[1],size[2]])
                coronal= dict.fromkeys(range(0, meas[i]),dimensions)
            if (i==1):
                dimensions=zeros([size[0],size[2]])
                sagital= dict.fromkeys(range(0, meas[i]),dimensions)
            if (i==2):
                dimensions=zeros([size[0],size[1]])
                axial= dict.fromkeys(range(0, meas[i]),dimensions)
        
        ccl= meas[0]+meas[1]+meas[2]
        
        im = {'corrcoef':zeros([1,ccl*10]), 'greatest':zeros([1,ccl]), 'template':zeros([1,ccl])}
        
        for i in range (0, meas[0]): #Gets the number of slides choosen
            a= beginning[0]+i #Starts in the slide selected as beginning and ends passing the one selected as end
            coronal[i]=np.uint8(((recon.get_fdata()[a,:,:])[:,:,0])/4)  #Get the slide of the reconstruction image.
        
        for i in range (0, meas[1]):    
            a= beginning[1]+i #Starts in the slide selected as beginning and ends passing the one selected as end
            sagital[i]=np.uint8(((recon.get_fdata()[:,a,:])[:,:,0])/4)  #Get the slide of the reconstruction image.
        
        for i in range (0, meas[2]):    
            a= beginning[2]+i #Starts in the slide selected as beginning and ends passing the one selected as end
            axial[i]=np.uint8(((recon.get_fdata()[:,:,a])[:,:,0])/4)  #Get the slide of the reconstruction image.
                    
        number=22 #Stablishes the base to number the templates
        
        def mean2(value):
            mean2value=np.sum(value)/np.size(value)
            return mean2value
        
        def corr2(R, T):
            R=R-mean2(R)
            T=T-mean2(T)
            corr=((R*T).sum())/(math.sqrt((R*R).sum()*(T*T).sum()))
            return corr
       
        for i in range (0,10): #Says that the process will repeat for the 10 templates
            number=number+1 #The first template will be 23
            volume='csf-aligned%d.nii.gz' %number #Construct the name of the template volume that will be loaded
            volume=nib.load(volume) #Load the template volume
            meascor=(meas[0]*i)
            meassag=(meas[1]*i)+(meas[0]*9)
            measax=(meas[2]*i)+(meas[0]*9)+(meas[1]*9)
            
            for j in range (0,ccl):
                if (j<meas[0]):
                    t=meascor+j #Define the position in which the results will be stored\
                    a=beginning[0]+j; #Select the slide that will be taken
                    slide=volume.get_fdata()[a,:,:] #Loads the slide
                    
                if (j>=meas[0] and j<(meas[0]+meas[1])):
                    t=meassag+j #Define the position in which the results will be stored\
                    a=beginning[0]+j-meas[0]; #Select the slide that will be taken
                    slide=volume.get_fdata()[:,a,:] #Loads the slide
                
                if (j>=(meas[0]+meas[1])):
                    t=measax+j #Define the position in which the results will be stored\
                    a=beginning[0]+j-meas[0]-meas[1]; #Select the slide that will be taken
                    slide=volume.get_fdata()[:,:,a] #Loads the slide    
                    
               #normalize the slide
                slide=np.uint8(256*(slide-slide.min())/(slide.max()-slide.min()))
               
                index=np.nonzero(slide) #get the positions of the non zero values 
                csfi=slide[np.nonzero(slide)] #get the values of the no zero values
              
                temp=csfi.shape #gets the number of non zero values 
                temp=temp[0]
                reconi=[]
                
                if (j<meas[0]):
                    for n in range (0,temp): #get the same indexes of the recon image
                        slide=coronal[j].item(index[0][n],index[1][n])
                        reconi.append(slide)
                                  
                if (j>=meas[0] and j<(meas[0]+meas[1])):
                    for n in range (0,temp): #get the same indexes of the recon image
                        jj=j-meas[0]
                        slide=sagital[jj].item(index[0][n],index[1][n])
                        reconi.append(slide)
                        
                if (j>=(meas[0]+meas[1])):
                    for n in range (0,temp): #get the same indexes of the recon image
                        jj=j-meas[0]-meas[1]
                        slide=axial[jj].item(index[0][n],index[1][n])
                        reconi.append(slide)
                        
                im['corrcoef'][0,t]=corr2(reconi, csfi)
             
            a=0
                 
        for a in range (0,ccl):
            im['greatest'][0,a]=-5
                
            for i in range (0,10):
                if (a<meas[0]):
                    j=(meas[0]*i)+a
                                    
                if (a>=meas[0] and a<(meas[1]+meas[0]) ):
                    j=(meas[1]*i)+(meas[0]*9)+a
                    
                if (a>=(meas[1]+meas[0])):    
                    j=(meas[2]*i)+(meas[0]*9)+(meas[1]*9)+a            
                    
                if im['corrcoef'][0,j]>im['greatest'][0,a]:
                    im['greatest'][0,a]=im['corrcoef'][0,j]
                    im['template'][0,a]=i+23
        
        t=im['template']
        t=t[0]
        a=t.tolist()
        j=im['greatest']
        jj=max(j[0])-((max(j[0]))/2.5)
        
        for i in range (0,ccl):
            if (j[0][i]>jj):
                a.append(im['template'][0,i])
        
        n=np.histogram(t,bins=[23,24,25,26,27,28,29,30,31,32,33])
        temp=np.argsort(n[0])[::-1]    
        tempi=temp+23
        temp=tempi[0]
        
        os.environ['temp']=str(temp)
       
        os.system('convert_xfm -omat InvAligned-$temp.xfm -inverse Temp-Recon-7dof-$temp.xfm')
        
        os.system('flirt -in recon.nii -ref ~/templates/template-$temp/template-$temp.nii \
        -out InvAligned-$temp.nii.gz -init InvAligned-$temp.xfm -applyxfm')
        
        os.chdir(path)
         
        os.system('cp alignment/InvAligned-$temp.nii.gz segmentation')
        os.system('cp alignment/InvAligned-$temp.xfm segmentation')
        
        os.chdir('segmentation')
        
        os.system('mv InvAligned-$temp.nii.gz Recon_final.nii.gz')
        os.system('mv InvAligned-$temp.xfm Recon_final.xfm')
    

    def __init__(self, *args, **kwargs):
        """
        The constructor
        """
        print("hello!!!")

        self.myName     = "Lizeth"
        self.inputDir   = ''
        self.ouptutDir  = ''

        pudb.set_trace()

        for k,v in kwargs.items():
            if k == 'inputDir':     self.inputDir   = v
            if k == 'outputDir':    self.outputDir  = v

    def run(self, options):
        """
        Define the code to be run by this plugin app.
        """
        print(Gstr_title)
        print('Version: %s' % self.get_version())
        os.chdir('/neuro/users/lizeth.machado/Public')
        #delete_ima()
        ima=self.verification_image()
        #brain_masking_tool()
        ima_mask=self.verification_image_mask()
        self.verification_image_mask_traspond()
        self.bias_field_correction(ima, ima_mask)
        self.writing_script(ima)
        self.alignment()


