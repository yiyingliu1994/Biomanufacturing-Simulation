
import SimFunctions
import SimClasses
import pandas
import numpy as np
import math
from sklearn.linear_model import LinearRegression
from scipy import stats
from matplotlib import pyplot as plt


np.random.seed(seed=123)


class regression():
	def __init__(self, FerData, Time):
		self.FerData = FerData
		self.Time = Time
		self.grate = 0
		self._P = 0
		self.sig_P = 0
		self.alpha = 0
		self._I = 0
		self.sig_I = 0

	def bootstrap(self):
		index = np.random.choice(len(self.FerData.iloc[:, 2]), len(self.FerData.iloc[:, 2]))
		self.FerData = self.FerData.iloc[index, :]

	def reg(self):
		y = np.log(self.FerData.iloc[:, 2]) - np.log(self.FerData.iloc[:, 1])
		y = y.values.reshape(-1, 1)
		X = np.array([self.Time for i in range(len(y))])
		X = X.reshape(-1, 1)
		clf = LinearRegression(fit_intercept=False)
		clf.fit(X, y)
		y_hat = clf.predict(X)
		self.grate = round(clf.coef_[0][0], 4)
		self.P = y - y_hat
		self.sig_P = round(np.std(self.P, ddof=1), 4)
		y2 = self.FerData.iloc[:, 3] / self.FerData.iloc[:, 2]
		y2 = y2.values.reshape(-1, 1)
		self.alpha = round(pow(np.prod(y2), 1.0 / len(y2)), 4)
		y_hat2 = self.FerData.iloc[:, 2] * self.alpha
		y_hat2 = y_hat2.values.reshape(-1, 1)
		exp_I = self.FerData.iloc[:, 3].values.reshape(-1, 1) / y_hat2
		self._I = np.log(exp_I)
		self.sig_I = round(np.std(self._I, ddof=1), 4)



class Batch(SimClasses.Entity):
    def __init__(self, name=""):
        super(Batch, self).__init__()
        self.Name = name
        self.Biomass = 0.0  # biomass for pre-culture
        self.External = False  # whether use external media (default = False)
        self.Protein = 0.0  # protein for antigen
        self.Impurity = 0.0  # impurities for antigen
        self.ChrCount = 0  # count number of chromatography steps


class Media(SimClasses.Entity):
    def __init__(self, name="", lifetime=2.0):
        super(Media, self).__init__()
        self.Name = name
        self.LifeTime = lifetime
        self.Process = "Ino"  # "Ino"/"Fer" for inoculum or main fermentation


def ECDF(Data):
    mass = Data.iloc[:, 1]
    mass = mass.tolist()
    mass.sort()
    U = np.random.uniform(0,1)
    m = len(mass)
    i = math.ceil((m - 1) * U)
    EmpDis = mass[i - 1] + (m - 1) * (mass[i] - mass[i - 1]) * (U - (i - 1) / (m - 1))
    return EmpDis


InoTime = {"a": 24.0, "b": 12.0, "clean": 1.5}
MedTime = {"a": 0.0, "b": 0.5, "clean": 0.6}
FerTime = {"a": 72.0, "b": 48.0, "clean": 6.0}
CenTime = {"a": 2.5, "b": 2.0, "clean": 0.2}
ChrTime = {"a": 8.0, "b": 7.0, "clean": 1.5}
FilTime = {"a": 2.0, "b": 2.0, "clean": 0.5}
QCTime = {"a": 2.0, "b": 2.0}


PreCulOrders = SimClasses.FIFOQueue()  # preculture buffer
InoPreCul_b = SimClasses.FIFOQueue()  # "b" preculture already start Inomedia prep (Priority to Orders queue) buffer
InoMedia = SimClasses.FIFOQueue()  # media use for inoculum
FerQueue = SimClasses.FIFOQueue()  # main fermenation
MainMedia = SimClasses.FIFOQueue()  # media use for main fermentation 
MediaToPrepare = SimClasses.FIFOQueue()  # record media preparation signals haven't been processed
DSPBuffer = SimClasses.FIFOQueue()  # downstream releasing


InoEqp = SimClasses.Resource()
MedEqp = SimClasses.Resource()
FerEqp = SimClasses.Resource()
CenEqp = SimClasses.Resource()
ChrEqp = SimClasses.Resource()
FilEqp = SimClasses.Resource()


InoEqp.SetUnits(5)
MedEqp.SetUnits(5)
FerEqp.SetUnits(5)
CenEqp.SetUnits(2)
ChrEqp.SetUnits(5)
FilEqp.SetUnits(2)


BusyUSP = SimClasses.CTStat()
BusyDSP = SimClasses.CTStat()


FerImp = 0.85  # impurity level at Fermentation
FinImp = 0.2  # threshold for QC

# Batch limit (we tried 3 different combinations: 4&2, 5&2, 6&3, 8&4)
USPLimit = 5
DSPLimit = 2


# Estimate Fermentation coefficients
AFerData = pandas.read_excel("Fer_A.xlsx")
BFerData = pandas.read_excel('Fer_B.xlsx')

A_reg = regression(AFerData, Time = FerTime['a'])
A_reg.reg()

B_reg = regression(BFerData, Time = FerTime['b'])
B_reg.reg()



# Centrifuge Distribution (uniform)
Cen_low = 0.4
Cen_high = 0.5


# Chromatography Distribution (uniform)
ChrData = pandas.read_excel('Chr.xlsx')

Qp_low = min(ChrData.iloc[:,1])
Qp_high = max(ChrData.iloc[:,1])

Qi_low = min(ChrData.iloc[:,2])
Qi_high = max(ChrData.iloc[:,2])

# Filtrition Distribution (uniform)
fil_low = 0.99
fil_high = 1.0



# record number of working antigen in USP Process, for USP releasing
NumUSP = 0

# record number of working antigen in Centrifuge & Chromatography Process, for DSP releasing
NumCenChr = 0

Calendar = SimClasses.EventCalendar()
batch_results = []
RunLength = 200
WarmUp = 500
Total = RunLength + WarmUp
Rep = 100  # number of replications

#Record input Biomass quality
biomass_input_a = []
biomass_input_b = []

def Arrival():
    global NumUSP  # number of batches in upstream (used for release policy)
    global Total
    global biomass_input_a
    global biomass_input_b
    cur_batch = Batch()  # generate a new entity
    Total = Total - 1

    if Total == RunLength:
        SimFunctions.Schedule(Calendar, "ClearIt", 0.000001)

    # schedule next arrival
    if Total > 0:
        SimFunctions.Schedule(Calendar, "Arrival", np.random.exponential(24.0))

    # type "a"
    if np.random.uniform(0.0, 1.0) < 0.25:
        # if it is type "a". Create the attributes for this entity a
        cur_batch.Name = "a"
        # use empirical distribution, convert to continuous
        cur_batch.Biomass = ECDF(AFerData)
        #cur_batch.Biomass = np.random.uniform(low = np.min(AFerData.iloc[:,1]), high = np.max(AFerData.iloc[:,1]))
        #cur_batch.Biomass = np.random.normal(loc= np.mean(AFerData.iloc[:,1]), scale=np.std(AFerData.iloc[:,1]), )
        biomass_input_a.append(cur_batch.Biomass)
        cur_batch.External = True

        # check the release policy; see Inoculum available or not
        if InoEqp.Busy == InoEqp.NumberOfUnits or NumUSP >= USPLimit:
            PreCulOrders.Add(cur_batch)
        else:
            NumUSP = NumUSP + 1
            BusyUSP.Record(float(NumUSP))
            InoEqp.Seize(1)
            SimFunctions.SchedulePlus(Calendar, "EndIno", InoTime[cur_batch.Name], cur_batch)
    else:
        # if it is type "b". Create the attributes for this entity, b
        cur_batch.Name = "b"
        # use empirical distribution, convert to continuous
        cur_batch.Biomass = ECDF(BFerData)
        #cur_batch.Biomass = np.random.uniform(low = np.min(BFerData.iloc[:,1]), high = np.max(BFerData.iloc[:,1]))
        #cur_batch.Biomass = np.random.normal(loc= np.mean(BFerData.iloc[:,1]), scale=np.std(BFerData.iloc[:,1]), )
        biomass_input_b.append(cur_batch.Biomass)
        cur_batch.External = False

        # see MediaPrep available or not
        if MedEqp.Busy == MedEqp.NumberOfUnits or NumUSP >= USPLimit:
            PreCulOrders.Add(cur_batch)
        else:
            NumUSP = NumUSP + 1
            BusyUSP.Record(float(NumUSP))
            MedEqp.Seize(1)
            InoPreCul_b.Add(cur_batch)  # put into the place to wait for media ready
            SimFunctions.SchedulePlus(Calendar, "EndInoMed", MedTime[cur_batch.Name], cur_batch.Name)



def EndIno(cur_batch):
    # end of inoculum trigger the cleaning of equipment
    SimFunctions.Schedule(Calendar, "EndCleanIno", InoTime["clean"])

    if cur_batch.External:
        # type "a" use external
        if FerEqp.Busy == FerEqp.NumberOfUnits:
            FerQueue.Add(cur_batch)
        else:
            FerEqp.Seize(1)
            SimFunctions.SchedulePlus(Calendar, "EndFer", FerTime[cur_batch.Name], cur_batch)
    else:
        # type "b"
        if FerEqp.Busy == FerEqp.NumberOfUnits:
            FerQueue.Add(cur_batch)
        else:
            while MainMedia.NumQueue() > 0:
                cur_media = MainMedia.Remove()  # remove the media from media preparation for main fermentation
                if SimClasses.Clock <= (cur_media.CreateTime + cur_media.LifeTime):  # not expired
                    # Suppose we have two separate queue: MainMedia, InoMedia
                    FerEqp.Seize(1)
                    SimFunctions.SchedulePlus(Calendar, "EndFer", FerTime[cur_batch.Name], cur_batch)
                    break
            # otherwise, all media for main fermentation expired/removed, or no available media, schedule
            if MedEqp.Busy == MedEqp.NumberOfUnits:
                media_signal = Media(name=cur_batch.Name)
                media_signal.Process = "Fer"
                MediaToPrepare.Add(media_signal)
            else:
                MedEqp.Seize(1)
                SimFunctions.SchedulePlus(Calendar, "EndMainMed", MedTime[cur_batch.Name], cur_batch.Name)




def EndInoMed(this_name):  # end of innoculum media preparation
    cur_media = Media(name=this_name)  # generate a new entity of media for innoclume requirement
    cur_media.Process = "Ino"
    SimFunctions.Schedule(Calendar, "EndCleanMed", MedTime["clean"])  # trigger cleaning

    if InoMedia.NumQueue() > 0:
        InoMedia.Add(cur_media)
    elif InoEqp.Busy == InoEqp.NumberOfUnits:
        InoMedia.Add(cur_media)  # For innoclum, we have three queue: (1)PreCulOrders;
        # (2) InoPreCul_b (it is along with InoMedia);
    else:
        InoEqp.Seize(1)
        cur_batch = InoPreCul_b.Remove()
        SimFunctions.SchedulePlus(Calendar, "EndIno", InoTime[cur_batch.Name], cur_batch)

        # trigger the media preparation for main fermentation
        MainMediaPrepTime = InoTime[cur_batch.Name] - MedTime[cur_batch.Name] - 0.00001
        SimFunctions.SchedulePlus(Calendar, "StartMainMed", MainMediaPrepTime, cur_batch.Name)




def EndCleanIno():
    global NumUSP

    if InoMedia.NumQueue() > 0:
        while InoMedia.NumQueue() > 0:
            cur_media = InoMedia.Remove()
            if SimClasses.Clock <= (cur_media.CreateTime + cur_media.LifeTime):  # not expired
                cur_batch = InoPreCul_b.Remove()
                SimFunctions.SchedulePlus(Calendar, "EndIno", InoTime[cur_batch.Name], cur_batch)

                # trigger the media preparation for main fermentation
                MainMediaPrepTime = InoTime[cur_batch.Name] - MedTime[cur_batch.Name] - 0.00001
                SimFunctions.SchedulePlus(Calendar, "StartMainMed", MainMediaPrepTime, cur_batch.Name)
                return
        # otherwise, all media for inoculum expired/removed, schedule
        cur_batch = InoPreCul_b.ThisQueue[0]
        if MedEqp.Busy == MedEqp.NumberOfUnits:
            media_signal = Media(name=cur_batch.Name)
            media_signal.Process = "Ino"
            MediaToPrepare.Add(media_signal)
        else:
            MedEqp.Seize(1)
            SimFunctions.SchedulePlus(Calendar, "EndInoMed", MedTime[cur_batch.Name], cur_batch.Name)
        # check orders list
        if PreCulOrders.NumQueue() > 0 and NumUSP < USPLimit:
            cur_batch = PreCulOrders.ThisQueue[0]
            if cur_batch.External:
                NumUSP = NumUSP + 1
                BusyUSP.Record(float(NumUSP))
                cur_batch = PreCulOrders.Remove()
                SimFunctions.SchedulePlus(Calendar, "EndIno", InoTime[cur_batch.Name], cur_batch)
            else:
                InoEqp.Free(1)
    elif PreCulOrders.NumQueue() > 0 and NumUSP < USPLimit:
        cur_batch = PreCulOrders.ThisQueue[0]
        if cur_batch.External:
            NumUSP = NumUSP + 1
            BusyUSP.Record(float(NumUSP))
            cur_batch = PreCulOrders.Remove()
            SimFunctions.SchedulePlus(Calendar, "EndIno", InoTime[cur_batch.Name], cur_batch)
        else:
            InoEqp.Free(1)
    else:
        InoEqp.Free(1)




def EndCleanMed():
    global NumUSP
    if MediaToPrepare.NumQueue() > 0:  # MediaToPrepare: the record of media to prepare
        media_signal = MediaToPrepare.Remove()
        if media_signal.Process == "Ino":  # prepare for InoMedia
            SimFunctions.SchedulePlus(Calendar, "EndInoMed", MedTime[media_signal.Name], media_signal.Name)
        else:  # prepare for MainMedia
            SimFunctions.SchedulePlus(Calendar, "EndMainMed", MedTime[media_signal.Name], media_signal.Name)
    elif PreCulOrders.NumQueue() > 0 and NumUSP < USPLimit:
        cur_batch = PreCulOrders.ThisQueue[0]
        if not cur_batch.External:  # if not using external media
            NumUSP = NumUSP + 1
            BusyUSP.Record(float(NumUSP))
            cur_batch = PreCulOrders.Remove()
            InoPreCul_b.Add(cur_batch)
            SimFunctions.SchedulePlus(Calendar, "EndInoMed", MedTime[cur_batch.Name], cur_batch.Name)
        else:
            MedEqp.Free(1)
    else:
        MedEqp.Free(1)




def StartMainMed(this_name):  # start the media preparation for main fermentation
    if MedEqp.Busy == MedEqp.NumberOfUnits:
        media_signal = Media(name=this_name)
        media_signal.Process = "Fer"
        MediaToPrepare.Add(media_signal)
    else:
        MedEqp.Seize(1)
        SimFunctions.SchedulePlus(Calendar, "EndMainMed", MedTime[this_name], this_name)




def EndMainMed(this_name):
    cur_media = Media(name=this_name)
    cur_media.Process = "Fer"
    SimFunctions.Schedule(Calendar, "EndCleanMed", MedTime["clean"])
    if MainMedia.NumQueue() > 0:
        MainMedia.Add(cur_media)
    elif FerEqp.Busy == FerEqp.NumberOfUnits:
        MainMedia.Add(cur_media)  # wait for main fermentation
    elif FerQueue.NumQueue() > 0:
        cur_batch = FerQueue.ThisQueue[0]
        if not cur_batch.External:  # if not using external media (actually must be true)
            cur_batch = FerQueue.Remove()
            FerEqp.Seize(1)
            SimFunctions.SchedulePlus(Calendar, "EndFer", FerTime[cur_batch.Name], cur_batch)
        else:
            MainMedia.Add(cur_media)  # not really executed
    else:
        MainMedia.Add(cur_media)  # not really executed




def EndFer(cur_batch):
    global batch_results
    global NumCenChr
    global NumUSP
    # USP release or not
    if PreCulOrders.NumQueue() > 0:
        next_batch = PreCulOrders.ThisQueue[0]
        if next_batch.External:
            # use external media
            if InoEqp.Busy < InoEqp.NumberOfUnits:
                next_batch = PreCulOrders.Remove()
                InoEqp.Seize(1)
                SimFunctions.SchedulePlus(Calendar, "EndIno", InoTime[next_batch.Name], next_batch)
            else:
                NumUSP = NumUSP - 1
                BusyUSP.Record(float(NumUSP))
        else:
            if MedEqp.Busy < MedEqp.NumberOfUnits:
                next_batch = PreCulOrders.Remove()
                MedEqp.Seize(1)
                InoPreCul_b.Add(next_batch)
                SimFunctions.SchedulePlus(Calendar, "EndInoMed", MedTime[next_batch.Name], next_batch.Name)
            else:
                NumUSP = NumUSP - 1
                BusyUSP.Record(float(NumUSP))
    else:
        NumUSP = NumUSP - 1
        BusyUSP.Record(float(NumUSP))

    # Fermentation output
    if cur_batch.Name == "a":
        cur_batch.Protein = cur_batch.Biomass * np.exp(A_reg.grate * FerTime['a'] + np.random.normal(0.0, A_reg.sig_P))
        cur_batch.Impurity = A_reg.alpha * cur_batch.Protein * np.exp(np.random.normal(0.0, A_reg.sig_I))
    else:
        cur_batch.Protein = cur_batch.Biomass * np.exp(B_reg.grate * FerTime['b'] + np.random.normal(0.0, B_reg.sig_P))
        cur_batch.Impurity = B_reg.alpha * cur_batch.Protein * np.exp(np.random.normal(0.0, B_reg.sig_I))

    # schedule clean
    SimFunctions.Schedule(Calendar, "EndCleanFer", FerTime["clean"])

    # compute the threshold
    cur_p = cur_batch.Protein
    cur_i = cur_batch.Impurity
    potential_imp = cur_i / (cur_i + cur_p)

    # check if can meet final quality or drop batch
    if potential_imp > FerImp:
        # Drop the batch
        cycletime = SimClasses.Clock - cur_batch.CreateTime
        batch_results.append([cur_batch.Name, cur_batch.Protein, cur_batch.Impurity, cycletime,
                              cur_batch.CreateTime, SimClasses.Clock, "FerDropped"])
    else:
        # see DSP availability
        if CenEqp.Busy == CenEqp.NumberOfUnits or NumCenChr >= DSPLimit:
            DSPBuffer.Add(cur_batch)
        else:
            CenEqp.Seize(1)
            NumCenChr = NumCenChr + 1
            BusyDSP.Record(float(NumCenChr))
            SimFunctions.SchedulePlus(Calendar, "EndCen", CenTime[cur_batch.Name], cur_batch)




def EndCleanFer():
    if FerQueue.NumQueue() > 0:
        cur_batch = FerQueue.ThisQueue[0]
        if cur_batch.External:  # if use external media
            cur_batch = FerQueue.Remove()
            SimFunctions.SchedulePlus(Calendar, "EndFer", FerTime[cur_batch.Name], cur_batch)
        else:  # if use internal media
            while MainMedia.NumQueue() > 0:
                cur_media = MainMedia.Remove()  # remove the media from media preparation for main fermentation
                if SimClasses.Clock <= (cur_media.CreateTime + cur_media.LifeTime):  # not expired
                    SimFunctions.SchedulePlus(Calendar, "EndFer", FerTime[cur_batch.Name], cur_batch)
                    return

            FerEqp.Free(1)
            # otherwise, all media for main fermentation expired/removed, or no available media, schedule
            if MedEqp.Busy == MedEqp.NumberOfUnits:
                media_signal = Media(name=cur_batch.Name)
                media_signal.Process = "Fer"
                MediaToPrepare.Add(media_signal)
            else:
                MedEqp.Seize(1)
                SimFunctions.SchedulePlus(Calendar, "EndMainMed", MedTime[cur_batch.Name], cur_batch.Name)
    else:
        FerEqp.Free(1)




def EndCen(cur_batch):
    global batch_results
    global NumCenChr
    Q = np.random.uniform(Cen_low, Cen_high)
    cur_batch.Impurity = Q * cur_batch.Impurity

    SimFunctions.Schedule(Calendar, "EndCleanCen", CenTime["clean"])

    prop_imp = cur_batch.Impurity / (cur_batch.Impurity + cur_batch.Protein)
    if prop_imp <= FinImp:  # skip chromatography, directly go to filtration
        if FilEqp.Busy == FilEqp.NumberOfUnits:  # violate no-wait constraint
            # drop the batch and record
            cycletime = SimClasses.Clock - cur_batch.CreateTime
            batch_results.append([cur_batch.Name, cur_batch.Protein, cur_batch.Impurity, cycletime,
                                  cur_batch.CreateTime, SimClasses.Clock, "DSPDropped"])
        else:
            FilEqp.Seize(1)
            SimFunctions.SchedulePlus(Calendar, "EndFil", FilTime[cur_batch.Name], cur_batch)
    else:
        if ChrEqp.Busy == ChrEqp.NumberOfUnits:  # violate no-wait constraint
            # drop the batch and record
            cycletime = SimClasses.Clock - cur_batch.CreateTime
            batch_results.append([cur_batch.Name, cur_batch.Protein, cur_batch.Impurity, cycletime,
                                  cur_batch.CreateTime, SimClasses.Clock, "DSPDropped"])
        else:
            ChrEqp.Seize(1)
            SimFunctions.SchedulePlus(Calendar, "EndChr", FilTime[cur_batch.Name], cur_batch)
    # check to release next batch
    if DSPBuffer.NumQueue() > 0 and CenEqp.Busy < CenEqp.NumberOfUnits:
        next_batch = DSPBuffer.Remove()
        CenEqp.Seize(1)
        SimFunctions.SchedulePlus(Calendar, "EndCen", CenTime[next_batch.Name], next_batch)
    else:
        NumCenChr = NumCenChr - 1
        BusyDSP.Record(float(NumCenChr))




def EndCleanCen():
    global NumCenChr
    if DSPBuffer.NumQueue() > 0 and NumCenChr < DSPLimit:
        cur_batch = DSPBuffer.Remove()
        NumCenChr = NumCenChr + 1
        BusyDSP.Record(float(NumCenChr))
        SimFunctions.SchedulePlus(Calendar, "EndCen", CenTime[cur_batch.Name], cur_batch)
    else:
        CenEqp.Free(1)


def EndChr(cur_batch):
	global NumCenChr
	step = cur_batch.ChrCount
	SimFunctions.Schedule(Calendar, "EndCleanChr", ChrTime["clean"])
	if step < 3: # change this 
		prop_imp = cur_batch.Impurity / (cur_batch.Impurity +cur_batch.Protein)
		if prop_imp <= FinImp:  # go to filtration
			if FilEqp.Busy == FilEqp.NumberOfUnits:  # violate no-wait constraint
				# drop the batch and record
				cycletime = SimClasses.Clock - cur_batch.CreateTime
				batch_results.append([cur_batch.Name, cur_batch.Protein, cur_batch.Impurity, cycletime,
									  cur_batch.CreateTime, SimClasses.Clock, "DSPDropped"])
			else:
				FilEqp.Seize(1)
				SimFunctions.SchedulePlus(Calendar, "EndFil", FilTime[cur_batch.Name], cur_batch)
		else:
			if ChrEqp.Busy == ChrEqp.NumberOfUnits:  # violate no-wait constraint
				# drop the batch and record
				cycletime = SimClasses.Clock - cur_batch.CreateTime
				batch_results.append([cur_batch.Name, cur_batch.Protein, cur_batch.Impurity, cycletime,
									  cur_batch.CreateTime, SimClasses.Clock, "DSPDropped"])
				# check to release next batch
			else:
				ChrEqp.Seize(1)
				# compute corresponding pooling window statistics
				Qp = np.random.uniform(Qp_low, Qp_high)
				Qi = np.random.uniform(Qi_low, Qi_high)
				cur_batch.Protein = Qp * cur_batch.Protein
				cur_batch.Impurity = Qi * cur_batch.Impurity
				cur_batch.ChrCount = cur_batch.ChrCount + 1  # chromatography time
				SimFunctions.SchedulePlus(Calendar, "EndChr", ChrTime[cur_batch.Name], cur_batch)

			if DSPBuffer.NumQueue() > 0 and CenEqp.Busy < CenEqp.NumberOfUnits:
				next_batch = DSPBuffer.Remove()
				CenEqp.Seize(1)
				SimFunctions.SchedulePlus(Calendar, "EndCen", CenTime[next_batch.Name], next_batch)
			else:
				NumCenChr = NumCenChr - 1
				BusyDSP.Record(float(NumCenChr))

		# check to release next batch
		if DSPBuffer.NumQueue() > 0 and CenEqp.Busy < CenEqp.NumberOfUnits:
			next_batch = DSPBuffer.Remove()
			CenEqp.Seize(1)
			SimFunctions.SchedulePlus(Calendar, "EndCen", CenTime[next_batch.Name], next_batch)
		else:
			NumCenChr = NumCenChr - 1
			BusyDSP.Record(float(NumCenChr))


	else:
	# finish chromatography for this batch
		if FilEqp.Busy == FilEqp.NumberOfUnits:  # violate no-wait constraint
			# drop the batch and record
			cycletime = SimClasses.Clock - cur_batch.CreateTime
			batch_results.append([cur_batch.Name, cur_batch.Protein, cur_batch.Impurity, cycletime,
								  cur_batch.CreateTime, SimClasses.Clock, "DSPDropped"])
		else:
			FilEqp.Seize(1)
			SimFunctions.SchedulePlus(Calendar, "EndFil", FilTime[cur_batch.Name], cur_batch)
		# check to release next batch
		if DSPBuffer.NumQueue() > 0 and CenEqp.Busy < CenEqp.NumberOfUnits:
			next_batch = DSPBuffer.Remove()
			CenEqp.Seize(1)
			SimFunctions.SchedulePlus(Calendar, "EndCen", CenTime[next_batch.Name], next_batch)
		else:
			NumCenChr = NumCenChr - 1
			BusyDSP.Record(float(NumCenChr))


def EndCleanChr():
    ChrEqp.Free(1)


def EndFil(cur_batch):
    Q = np.random.uniform(fil_low, fil_high)
    cur_batch.Impurity = Q * cur_batch.Impurity
    SimFunctions.Schedule(Calendar, "EndCleanFil", FilTime["clean"])
    SimFunctions.SchedulePlus(Calendar, "EndQC", QCTime[cur_batch.Name], cur_batch)


def EndCleanFil():
    FilEqp.Free(1)



def EndQC(cur_batch):
    global batch_results
    cycletime = SimClasses.Clock - cur_batch.CreateTime
    if cur_batch.Impurity / (cur_batch.Impurity + cur_batch.Protein) > FinImp:
        # Not qualified
        batch_results.append([cur_batch.Name, cur_batch.Protein, cur_batch.Impurity, cycletime,
                              cur_batch.CreateTime, SimClasses.Clock, "Not qualified"])
    else:
        batch_results.append([cur_batch.Name, cur_batch.Protein, cur_batch.Impurity, cycletime,
                              cur_batch.CreateTime, SimClasses.Clock, "Finished"])




TheCTStats = [BusyUSP, BusyDSP, InoEqp.NumBusy, FerEqp.NumBusy, CenEqp.NumBusy, ChrEqp.NumBusy, FilEqp.NumBusy]
TheDTStats = []
TheQueues = [PreCulOrders, InoPreCul_b, FerQueue, InoMedia, MainMedia, MediaToPrepare, DSPBuffer]
TheResources = [InoEqp, MedEqp, FerEqp, CenEqp, ChrEqp, FilEqp]

output_a = []
output_b = []
output_cycle = []


for reps in range(Rep):
    batch_results = []
    NumUSP = 0
    NumCenChr = 0
    Total = WarmUp + RunLength
    SimFunctions.SimFunctionsInit(Calendar, TheQueues, TheCTStats, TheDTStats, TheResources)
    SimFunctions.Schedule(Calendar, "Arrival", np.random.exponential(24.0))

    while Calendar.N() > 0:
        NextEvent = Calendar.Remove()
        SimClasses.Clock = NextEvent.EventTime
        if NextEvent.EventType == "Arrival":
            Arrival()
        elif NextEvent.EventType == "EndIno":
            EndIno(NextEvent.WhichObject)
        elif NextEvent.EventType == "EndInoMed":
            EndInoMed(NextEvent.WhichObject)
        elif NextEvent.EventType == "EndCleanIno":
            EndCleanIno()
        elif NextEvent.EventType == "EndCleanMed":
            EndCleanMed()
        elif NextEvent.EventType == "StartMainMed":
            StartMainMed(NextEvent.WhichObject)
        elif NextEvent.EventType == "EndMainMed":
            EndMainMed(NextEvent.WhichObject)
        elif NextEvent.EventType == "EndFer":
            EndFer(NextEvent.WhichObject)
        elif NextEvent.EventType == "EndCleanFer":
            EndCleanFer()
        elif NextEvent.EventType == "EndCen":
            EndCen(NextEvent.WhichObject)
        elif NextEvent.EventType == "EndCleanCen":
            EndCleanCen()
        elif NextEvent.EventType == "EndChr":
            EndChr(NextEvent.WhichObject)
        elif NextEvent.EventType == "EndCleanChr":
            EndCleanChr()
        elif NextEvent.EventType == "EndFil":
            EndFil(NextEvent.WhichObject)
        elif NextEvent.EventType == "EndCleanFil":
            EndCleanFil()
        elif NextEvent.EventType == "EndQC":
            EndQC(NextEvent.WhichObject)
        elif NextEvent.EventType == "ClearIt":
            SimFunctions.ClearStats(TheCTStats, TheDTStats)
        
    print ("This is No.%i" %reps, "replication")
    # print(batch_results)
    rep_output = np.asarray(batch_results)
    col = ["Name", "Protein", "Impurity", "CycleTime", "StartTime", "EndTime", "Status"]
    rep_output = pandas.DataFrame(data=rep_output, columns=col)
    rep_output[["StartTime"]] = rep_output[["StartTime"]].apply(pandas.to_numeric)
    rep_output = rep_output.sort_values("StartTime")
    rep_output = rep_output.iloc[WarmUp:, :]
    # compute statistics for each replication
    df_a = rep_output.loc[rep_output['Name'] == "a"]
    df_b = rep_output.loc[rep_output['Name'] == "b"]
    # a
    Finished_a = sum(list(df_a['Status'] == "Finished"))
    NotQua_a = sum(list(df_a['Status'] == "Not qualified"))
    FerDrop_a = sum(list(df_a['Status'] == "FerDropped"))
    DSPDrop_a = sum(list(df_a['Status'] == "DSPDropped"))
    df_a = df_a.loc[df_a['Status'] != "FerDropped"]
    df_a = df_a.loc[df_a['Status'] != "DSPDropped"]
    df_a = df_a.loc[df_a['Status'] != "Not qualified"]
    ProMean_a = df_a['Protein'].astype(float).mean()
    ImpMean_a = df_a['Impurity'].astype(float).mean()
    CycMean_a = df_a['CycleTime'].astype(float).mean()
    ProVar_a = df_a['Protein'].astype(float).var()
    ImpVar_a = df_a['Impurity'].astype(float).var()
    CycVar_a = df_a['CycleTime'].astype(float).var()
    ThrPut_a = ProMean_a / CycMean_a
    ProTotal_a = df_a['Protein'].astype(float).sum()
    Pass_a = Finished_a / (Finished_a + NotQua_a + FerDrop_a + DSPDrop_a) if (Finished_a + NotQua_a + FerDrop_a + DSPDrop_a) != 0 else 0
    # b
    Finished_b = sum(list(df_b['Status'] == "Finished"))
    NotQua_b = sum(list(df_b['Status'] == "Not qualified"))
    FerDrop_b = sum(list(df_b['Status'] == "FerDropped"))
    DSPDrop_b = sum(list(df_b['Status'] == "DSPDropped"))
    df_b = df_b.loc[df_b['Status'] != "FerDropped"]
    df_b = df_b.loc[df_b['Status'] != "DSPDropped"]
    df_b = df_b.loc[df_b['Status'] != "Not qualified"]
    ProMean_b = df_b['Protein'].astype(float).mean()
    ImpMean_b = df_b['Impurity'].astype(float).mean()
    CycMean_b = df_b['CycleTime'].astype(float).mean()
    ProVar_b = df_b['Protein'].astype(float).var()
    ImpVar_b = df_b['Impurity'].astype(float).var()
    CycVar_b = df_b['CycleTime'].astype(float).var()
    ThrPut_b = ProMean_b / CycMean_b
    ProTotal_b = df_b['Protein'].astype(float).sum()
    Pass_b = Finished_b / (Finished_b + NotQua_b + FerDrop_b + DSPDrop_b) if (Finished_b + NotQua_b + FerDrop_b + DSPDrop_b) != 0 else 0
    # utilization
    Ino_Utl = InoEqp.Mean() / InoEqp.NumberOfUnits
    Fer_Utl = FerEqp.Mean() / FerEqp.NumberOfUnits
    Cen_Utl = CenEqp.Mean() / CenEqp.NumberOfUnits
    Chr_Utl = ChrEqp.Mean() / ChrEqp.NumberOfUnits
    Fil_Utl = FilEqp.Mean() / FilEqp.NumberOfUnits
    output_a.append([ProMean_a, ImpMean_a, CycMean_a, ProVar_a, ImpVar_a, CycVar_a, ThrPut_a,
                     ProTotal_a, Pass_a, Ino_Utl, Fer_Utl, Cen_Utl, Chr_Utl, Fil_Utl, Finished_a,
                     NotQua_a, FerDrop_a, DSPDrop_a])
    output_b.append([ProMean_b, ImpMean_b, CycMean_b, ProVar_b, ImpVar_b, CycVar_b, ThrPut_b,
                     ProTotal_b, Pass_b, Ino_Utl, Fer_Utl, Cen_Utl, Chr_Utl, Fil_Utl, Finished_b,
                     NotQua_b, FerDrop_b, DSPDrop_b])
    


# summary output
col = ["ProteinMean", "ImpurityMean", "CycleTimeMean", "ProteinVar", "ImpurityVar", "CycleTimeVar", "Throughput",
       "TotalYield", "PassRate", "Ino_Utilization", "Fer_Utilization", "Cen_Utilization", "Chr_Utilization",
       "Fil_Utilization", "Finished #", "Not Qualified #", "FerDropped #", "DSPDropped #"]
output_a = np.asarray(output_a)
output_a = pandas.DataFrame(data=output_a, columns=col)

output_b = np.asarray(output_b)
output_b = pandas.DataFrame(data=output_b, columns=col)

output_a.to_csv("output_a_Chr_step_4.csv", sep=",")
output_b.to_csv("output_b_Chr_step_4.csv", sep=",")





