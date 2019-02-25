import mrcfile

def save_vol_as_map(fname, vol, vsizes, label):
    mrc = mrcfile.new(fname, overwrite=True, data=vol)
    mrc.header['cella']['x'] = np.array(vsizes[0])
    mrc.header['cella']['y'] = np.array(vsizes[1])
    mrc.header['cella']['z'] = np.array(vsizes[2])
    for i in range(ceil(len(label)/80.)):
        mrc.header['label'][i] = label[i*80:(i+1)*80]
    mrc.close()
