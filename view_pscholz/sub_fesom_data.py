# Patrick Scholz, 14.12.2017
import numpy as np
import time
from netCDF4 import Dataset
from set_inputarray import *
from sub_fesom_mesh import fesom_vector_rot
global inputarray
	
#+_____________________________________________________________________________+
#|                                                                             |
#|                        *** FESOM DATA CLASS ***                             |
#|                                                                             |
#+_____________________________________________________________________________+
class fesom_data(object):
	#____data run variables______________________
	var                         = ''
	runid, path, descript       = '', '', ''
	
	#____data time variables_____________________
	year,  month, record, depth = [], [], [], []
	str_time, str_dep           = '', ''
	
	#____data projection variables_______________
	proj, proj_lon, proj_lat    = '', 0.0, 90.0
	
	#____data description info___________________
	sname, lname, unit	        = '', '', ''
	
	#____data plot varaibles_____________________
	cmap, crange, cnumb         = 'grads', [], []
	which_plot                  = ''
	
	#____data varaible___________________________
	value, value2               = [], []
	which_mean                  = 'monthly'
	anom                        = False
	
	
	#___INIT DATA OBJECT_______________________________________________________#
	#
	#__________________________________________________________________________#
	def __init__(self,inputarray):
		if len(inputarray)!=0:
			self.runid		= inputarray['data_id']
			self.path		= inputarray['data_dir1']
			self.proj		= inputarray['proj']
			self.proj_lon	= inputarray['proj_lon']
			self.proj_lat	= inputarray['proj_lat']
			self.which_plot	= inputarray['which_plot']
	
	
#___LOAD FESOM2.0 DATA__________________________________________________________
#
#_______________________________________________________________________________
def fesom_load_data_horiz(mesh,data):
	
	#___________________________________________________________________________
	# number of years to average 
	nyi = data.year[1]-data.year[0]+1
	
	# array with years
	ayi = np.linspace(data.year[0],data.year[1],nyi,dtype='uint16')
	# number of month to average
	nmi = np.array(data.month).size
	# array with month
	ami = np.array(data.month)
	
	print('')
	print('     -----+-----------------------------------+------------')
	print('     Year |               MON                 |')
	print('     -----+-----------------------------------+------------')
	print('     --> '+data.path)
	print('     --> '+data.var)
	#____START YEAR LOOP________________________________________________________
	aux_datavar = data.var
	for yi in range(0, nyi):
		print('     {:4.0f} |'.format(ayi[yi])),
		tstart = time.time();
		if data.var.find('norm')==0 or data.var.find('vec')==0:
			aux_datavar = 'u'
			aux_datavar2 = 'v'
			#_______________________________________________________________________
			# build filename where data are located
			fname  =   data.path+'/'\
						+ aux_datavar + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
			fname2 =   data.path+'/'\
						+ aux_datavar2 + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
			#_______________________________________________________________________
			# open netcdf file --> read NETCDF4
			ncid  = Dataset(fname, 'r') 
			ncid2 = Dataset(fname2, 'r') 
		else:
			#_______________________________________________________________________
			# build filename where data are located
			fname =   data.path+'/'\
					+ data.var + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
			
			#_______________________________________________________________________
			# open netcdf file --> read NETCDF4
			ncid = Dataset(fname, 'r') 
		
		
		# read dimension of variable 
		ncdims=ncid.variables[aux_datavar].shape
		#print(ncdims)
		
		# number of dimension is 2d or 3d
		dim_num = len(ncdims)
		
		# check time dimension of loaded data aka number of records
		nrec  = ncdims[0]
		nsmple= ncdims[1]
		
		#_______________________________________________________________________
		# at first loaded year allocate matrices to calc means
		if yi==0:
			
			# if data are 3d and no loadable depth layers are given load all fesom 
			# layers
			if dim_num==3 and len(data.depth)==0:
				# allocate array for mean data
				mean_data = np.zeros((nsmple,ncdims[2]),dtype='float32')
			else:
				# allocate array for mean data
				mean_data = np.zeros((nsmple,),dtype='float32')
			
			if data.var.find('vec')==0:
				mean_data2 = np.zeros((nsmple),dtype='float32')
		
		#_______________________________________________________________________
		ncval = np.array(ncid.variables[aux_datavar][:])
		if data.var.find('norm')==0:
			ncval = np.sqrt(ncval[:]**2+np.array(ncid2.variables[aux_datavar2][:]**2))
		elif data.var.find('vec')==0:
			ncval2 = ncid2.variables[aux_datavar2][:]
		#_______________________________________________________________________
		# check if a single time slice record is selcted
		# no ?  --> than use data.month to select time subset
		if len(data.record)==0:
			#___________________________________________________________________
			if nrec==1:
				print('annual |'),
			elif nrec==12:
				print('monthly|'),
			#elif nrec==73:
				#print('5daily |'),
			elif nrec==365:
				print('daily  |'),
			else:	
				print(' --> error: temporal length of data unclear or not supported!!!|'),
				break
			#___________________________________________________________________
			# select monthly or seasonal subset of data
			if nmi<12 and nrec>=12:
				for ii in data.month: print('{:d}|'.format(ii)),
				#_______________________________________________________________
				# select from monthly data
				if nrec==12 : 
					ncval=fesom_time_depth_mean(mesh,ncval,np.array(data.month)-1,data.depth)
					if data.var.find('vec')==0: 
						ncval2=fesom_time_depth_mean(mesh,ncval2,np.array(data.month)-1,data.depth)
				#___________________________________________________________________
				# select from daily data
				elif nrec==365 : 
					# find out which day belongs to selected subset
					idx_sel = sel_timesubset_daily(data)
					ncval=fesom_time_depth_mean(mesh,ncval,idx_sel==True,data.depth)
					if data.var.find('vec')==0: 
						ncval2=fesom_time_depth_mean(mesh,ncval2,idx_sel==True,data.depth)
					del idx_sel
					
			#_______________________________________________________________________
			# select full data
			else:	
				if nrec>=12: 
					for ii in data.month: print('{:d}|'.format(ii)),
				ncval=fesom_time_depth_mean(mesh,ncval,[],data.depth)
				if data.var.find('vec')==0: 
					ncval2=fesom_time_depth_mean(mesh,ncval2,[],data.depth)
		#_______________________________________________________________________
		# yes ?  --> use data.record to select time slice
		else:
			if np.max(data.record)-1>nrec:
				print(' --> error: this record number doesnt exist in that file')
				break
			print('select single record time slice|'),
			for ii in data.record: print('{:d}|'.format(ii)),
			ncval=fesom_time_depth_mean(mesh,ncval,data.record[0]-1,data.depth)
			if data.var.find('vec')==0: 
				ncval2=fesom_time_depth_mean(mesh,ncval2,data.record[0]-1,data.depth)
		#_______________________________________________________________________
		# setup variable name, description and unit
		ncvar = ncid.variables[aux_datavar]
		#data.sname=ncvar.name
		data.sname=data.var
		for attrname in ncvar.ncattrs():
			if attrname=='description':
				data.lname=ncvar.getncattr(attrname)
			elif attrname=='units':
				data.unit='['+ncvar.getncattr(attrname)+']'
		#_______________________________________________________________________
		# close netcdf file
		ncid.close()
		
		#_______________________________________________________________________
		# now average over yearly subset
		mean_data = mean_data + ncval/nyi
		del ncval 
		if data.var.find('vec')==0: 
			mean_data2 = mean_data2 + ncval2/nyi
			del ncval2 
		
		#_______________________________________________________________________
		tend = time.time();		
		print(' --> t={:2.2f}s'.format(tend-tstart))
	#_______________________________________________________________________
	# in case of ice data set no ice to nan
	if (data.var=='a_ice'):
		mean_data=mean_data*100.0
	#if (data.var=='a_ice' or data.var=='m_ice'):
		#mean_data[mean_data==0.0]=np.nan
	
	#____ROTATE VECTORES FROM ROT TO GEO________________________________________
	# in future check if rotation is s till neccessary -> we want to/should write out 
	# the data already rotated ?!
	if data.var.find('vec')==0:
		mean_data,mean_data2=fesom_vector_rot(mesh,mean_data,mean_data2)
		
	#____END YEAR LOOP__________________________________________________________
	if nsmple==mesh.n2dn:
		# augment node data
		data.value = np.concatenate((mean_data,mean_data[mesh.pbndn_2d_i]))
		if data.var.find('vec')==0: 
			data.value2 = np.concatenate((mean_data2,mean_data2[mesh.pbndn_2d_i]))
	else:
		# augment elem data
		data.value = np.concatenate((mean_data,mean_data[mesh.pbndtri_2d_i]))
		if data.var.find('vec')==0: 
			data.value2 = np.concatenate((mean_data2,mean_data2[mesh.pbndtri_2d_i]))
			
		# delete trinangle that are seen as abnormal 
		if len(mesh.abnormtri_2d_i)!=0:
			data.value = np.delete(data.value,mesh.abnormtri_2d_i)
			if data.var.find('vec')==0: 
				data.value2 = np.delete(data.value2,mesh.abnormtri_2d_i)
	del mean_data
	
	#___________________________________________________________________________
	# set up descriptiv variable names
	if len(data.record)==0:
		mon_list   = np.array(['J','F','M','A','M','J','J','A','S','O','N','D'])
		mon_list2  = np.array(['January','February','March','April','May','June','July','August','September','October','Novemver','December'])
		# year info 
		if nyi == 1:
			str_time_1 = 'y: '+str(data.year[0])
		else:	
			str_time_1 = 'y: '+str(data.year[0])+'-'+str(data.year[1])
			
		# month info 
		if nmi == 1 and nrec>=12:
			str_time_1 = str_time_1+' ['+mon_list2[np.array(data.month)-1][0]+']'
		elif nmi>1 and nmi<12 and nrec>=12:
			str_mon=''
			for ii in range(0,nmi):
				if ami[ii]==12 : str_mon=mon_list[ami[ii]-1]+str_mon
				else 		   : str_mon=str_mon+mon_list[ami[ii]-1]
			str_time_1 = str_time_1 + ' [' + str_mon + ']'
		data.str_time = str_time_1
	else: 
		data.str_time = 'rec_i: '+str(data.record[0])
	
	if dim_num==3 and len(data.depth)!=0:
		if len(data.depth)<=1:
			data.str_dep=', dep: '+str(data.depth[0])+'m'
		else:
			data.str_dep=', dep: '+str(data.depth[0])+'m-'+str(data.depth[-1])+'m'
	#___________________________________________________________________________
	return data
	
	
#___DO VERTICAL INTERPOLATE AVERAGE OVER CERTAIN LAYERS_________________________
#
#_______________________________________________________________________________
def fesom_vinterp(data_in,mesh,levels):
	
	#___________________________________________________________________________
	# do vertical linear interpolation of certain layer
	dims = data_in.shape
	nsmple = dims[0] # is data_in defined on nodes or elements
	
	data_out = np.zeros((nsmple,),dtype='float32')
	aux_div  = np.zeros((nsmple,),dtype='float32')
	#for di in range(0,len(data.depth)):
	for di in range(0,len(levels)):
		# find upper and lower layer indices
		#idx_dwn = np.array(np.where( data.depth[di]<=abs(mesh.zlev))).squeeze()
		idx_dwn = np.array(np.where( levels[di]<=abs(mesh.zlev))).squeeze()
		idx_dwn = idx_dwn[0]
		idx_up  = idx_dwn-1
		if idx_up<0: idx_up=0
		
		# linear vertical interpolant
		deltaz   = abs(mesh.zlev[idx_dwn])-abs(mesh.zlev[idx_up])
		#deltaz_i = abs(mesh.zlev[idx_dwn])-data.depth[di]
		deltaz_i = abs(mesh.zlev[idx_dwn])-levels[di]
		
		# data_in is defined on nodes temp, salt, ssh ...
		if nsmple==mesh.n2dn:
			if deltaz_i==0:
				data_out = data_out + data_in[:,idx_dwn]
				aux_div[mesh.nodes_2d_iz[0:mesh.n2dn]>idx_dwn]=aux_div[mesh.nodes_2d_iz[0:mesh.n2dn]>idx_dwn]+1.0
				
			else:
				auxval     = data_in[:,idx_dwn]-(data_in[:,idx_dwn]-data_in[:,idx_up])*deltaz_i/deltaz
				auxval[mesh.nodes_2d_iz[0:mesh.n2dn]<idx_dwn]=0
				aux_div[mesh.nodes_2d_iz[0:mesh.n2dn]>idx_dwn]=aux_div[mesh.nodes_2d_iz[0:mesh.n2dn]>idx_dwn]+1.0
				data_out = data_out + auxval
				
		#data_in is defined on elements u,v....
		elif nsmple==mesh.n2de:
			if deltaz_i==0:
				data_out = data_out + data_in[:,idx_dwn]
				aux_div[mesh.elem0_2d_iz>idx_dwn]=aux_div[mesh.elem0_2d_iz>idx_dwn]+1.0
				
			else:
				auxval     = data_in[:,idx_dwn]-(data_in[:,idx_dwn]-data_in[:,idx_up])*deltaz_i/deltaz
				auxval[mesh.elem0_2d_iz<idx_dwn]=0
				aux_div[mesh.elem0_2d_iz>idx_dwn]=aux_div[mesh.elem0_2d_iz>idx_dwn]+1.0
				data_out = data_out + auxval
			
	#___________________________________________________________________________
	# do mean over averaged layers
	data_out = data_out/aux_div
	data_out[aux_div==0        ]=np.nan
	data_out[np.isinf(data_out)]=np.nan
	#___________________________________________________________________________
	return(data_out)
	
	
#___SELECT TIME SUBSET OF DAILY DATA TO CALCULATE SEASONAL MEAN_________________
#
#_______________________________________________________________________________
def sel_timesubset_daily(data):
	
	daypermonth=np.array([0,31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
	daypermonth=daypermonth.cumsum()
	idx_day2mon = np.zeros((365,))
	idx_day     = np.arange(1,365+1,1)
	for ii in range(0,12): 
		idx_day2mon[(idx_day>daypermonth[ii]) & (idx_day<=daypermonth[ii+1])]=ii+1
	del idx_day
	del daypermonth
				
	# find out which day belongs to selected subset
	idx_sel = np.in1d(idx_day2mon,data.month)
	del idx_day2mon	
	
	return(idx_sel)
	
	
#___SELECT SUBSET AND DO TEMPORAL AND VERTICAL AVERAGE__________________________
#
#_______________________________________________________________________________
def fesom_time_depth_mean(mesh,ncval,sel_time,levels):
	
	dim_num = len(ncval.shape)
	#___________________________________________________________________________
	# do no time selection
	if len(sel_time)==0:
		if dim_num==2:
			# case 2d variable
			ncval= ncval[:,:].mean(axis=0)
		elif dim_num==3:
			# case 3d variable
			ncval= ncval[:,:,:].mean(axis=0)
			# do vertical interpolation of certain layer if data.depth is not empty
			#if len(data.depth)!=0 : ncval=fesom_vinterp(ncval,mesh)
			if len(levels)!=0 : ncval=fesom_vinterp(ncval,mesh,levels)
		else:
			print(' --> error: more than 3 variable dimensions are not supported' )
	#___________________________________________________________________________
	# select certain time slice and average
	else:
		if dim_num==2:
			# case 3d variable
			ncval= ncval[sel_time,:].mean(axis=0)
		elif dim_num==3:
			# case 3d variable
			ncval= ncval[sel_time,:,:].mean(axis=0)
			# do vertical interpolation of certain layer if data.depth is not empty
			#if len(data.depth)!=0 : ncval=fesom_vinterp(ncval,mesh)
			if len(levels)!=0 : ncval=fesom_vinterp(ncval,mesh,levels)
		else:
			print(' --> error: more than 3 variable dimensions are not supported' )
	return(ncval)
	
	
#___CALCULATE ANOMALY BETWEEN TWO DATA OBJECTS__________________________________
#
#_______________________________________________________________________________
def fesom_data_anom(data,data2):
	anom = fesom_data([])
	
	#____data run variables______________________
	anom.var                            = data.var
	anom.runid, anom.path               = data.runid, data.path
	anom.descript                       = data2.descript+'-'+data.descript
	
	#____data time variables_____________________
	anom.year,  anom.month, anom.record = data.year, data.month, data.record
	anom.depth                          = data.depth
	if data.str_time!=data2.str_time:
		anom.str_time = data2.str_time+'-'+data.str_time
	else:
		anom.str_time = data.str_time
	anom.str_dep                        = data.str_dep
	
	#____data projection variables_______________
	anom.proj                           = data.proj
	anom.proj_lon, anom.proj_lat        = data.proj_lon, data.proj_lat
	
	#____data description info___________________
	anom.sname, anom.lname, anom.unit   = data.sname, data.lname, data.unit
	
	#____data plot varaibles_____________________
	#anom.cmap, anom.crange, anom.cnumb  = data.cmap, data.crange, data.cnumb
	anom.which_plot                     = data.which_plot 
	#if len(data.crange)!=0 :anom.crange=[]
	#____data varaible___________________________
	anom.value                          = data2.value-data.value
	anom.anom                           = True
	anom.value2                         = data.value2
	anom.which_mean                     = data.which_mean
	
	#___________________________________________________________________________
	return(anom)