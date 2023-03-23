# encoding:utf-8
"""
    @autor:yuexuezhi
    @contact:846924329@qq.com
    @file:
    @time:
    @decs:
"""

from osgeo import gdal,ogr
import time
import numpy as np
import os

try:
    progress = gdal.TermProgress_nocb
except:
    progress = gdal.TermProgress

def ext(in_layer,meta1):
    # 计算逆放射变换系数
    rev_tranform=gdal.InvGeoTransform(meta1[3])

    # datasource=ogr.Open(shp_dir)
    # in_layer=datasource.GetLayer(0)
    extent=in_layer.GetExtent()

    ulx,uly=map(round,gdal.ApplyGeoTransform(rev_tranform,extent[0],extent[3]))
    drx,dry=map(round,gdal.ApplyGeoTransform(rev_tranform,extent[1],extent[2]))

    # 求取重叠区域
    if ulx>meta1[0] or uly > meta1[1] or drx <= 0 or dry <= 0:
        print("there is no overlapping")
    else:
        # 列
        str_column =min(max(0,ulx),meta1[0]-1)
        end_column=max(min(drx,meta1[0]),0)
        # 行
        str_line =min(max(0,uly),meta1[1]-1)
        end_line=max(min(dry,meta1[1]),0)

        overlay_ext=[str_column,end_column,str_line,end_line]

    return overlay_ext

def filldata(nodata,out,code,layer,meta1,temfold,star_end,fidname):
    shpdriver = ogr.GetDriverByName("ESRI Shapefile")
    ref = layer.GetSpatialRef()
    # 创建字段
    # defn = layer.GetLayerDefn()
    # fieldIndex = defn.GetFieldIndex('percent')
    # if fieldIndex < 0:
    #     # 添加字段
    #     fieldDefn = ogr.FieldDefn('percent', ogr.OFSTFloat32)
    #     fieldDefn.SetPrecision(5)
    #     layer.CreateField(fieldDefn, 1)
    # fieldIndex2 = defn.GetFieldIndex('percent')
    #
    # if fieldIndex2 > 0:
    #     print("字段创建成功：", fieldIndex)

    field_name = ogr.FieldDefn(fidname, ogr.OFTReal)  ## 设置属性
    field_name.SetWidth(5)
    field_name.SetPrecision(5)## 设置长度
    layer.CreateField(field_name)  ## 创建字段

    for i in range(star_end[0],star_end[1],1):

        feature=layer.GetFeature(i)
        mem_dri = ogr.GetDriverByName('Memory')
        mem_ds = mem_dri.CreateDataSource(' ')
        out_layer = mem_ds.CreateLayer(' ', geom_type=ogr.wkbPolygon, srs=ref)
        out_layer.CreateFeature(feature)
        overlay_ext = ext(out_layer,meta1)

        ulx1, uly1 = gdal.ApplyGeoTransform(meta1[3], overlay_ext[0], overlay_ext[2])

        out_transform1 = (ulx1, meta1[3][1], 0, uly1, 0, meta1[3][5])

        xsize1 = overlay_ext[1] - overlay_ext[0]
        ysize1= overlay_ext[3] - overlay_ext[2]
        metadata = [xsize1, ysize1, meta1[2], out_transform1]
        ra_feat = vector2raster(out_layer, metadata)
        ra_feat_data = ra_feat.ReadAsArray()
        number_ra = np.sum(ra_feat_data)
        outmask=np.zeros((meta1[1],meta1[0]),np.uint8)
        outmask[overlay_ext[2]:overlay_ext[3],overlay_ext[0]:overlay_ext[1]]=ra_feat_data
    #
        mask_data=meta1[4][overlay_ext[2]:overlay_ext[3],overlay_ext[0]:overlay_ext[1]]
        mask = mask_data[np.where(ra_feat_data == 1)]
        number_leibie = np.sum(mask == code)

        #feature.SetField('20220310', number_leibie / number_ra)
        feature.SetField(fidname, 1-(number_leibie / number_ra))
        layer.SetFeature(feature)

    #     if number_leibie / number_ra > propotation:
    #         out[overlay_ext[2]:overlay_ext[3],overlay_ext[0]:overlay_ext[1]][np.where(ra_feat_data == 1)] = code
    #     else:
    #         out[overlay_ext[2]:overlay_ext[3],overlay_ext[0]:overlay_ext[1]][np.where(ra_feat_data == 1)] = nodata
    # return out

def vector2raster(layer1,metadata):
    """
矢量栅格化
    :param layer: 待栅格化的矢量图层
    :param metadata: 参考栅格的源数据
    :return: 与参考栅格一样大小的二值化图像
    """
    driver = gdal.GetDriverByName("MEM")
    outdataset = driver.Create('', metadata[0], metadata[1], 1, gdal.GDT_Byte)
    outdataset.SetGeoTransform(metadata[3])
    outdataset.SetProjection(metadata[2])
    gdal.RasterizeLayer(outdataset, [1], layer1, burn_values=[1])
    # RasterizeLayer 参数设置具体含义
    return outdataset

def calc_confidence_interval(msavi, ratio):

    bins = np.arange(start=msavi.min(), stop=msavi.max(), step=1)
    hist, bin_edges = np.histogram(msavi, bins=bins)

    nozero_index = np.where(bin_edges[:-1] != 0)
    zero_count = np.where(msavi == 0)[0].shape[0]
    hist = hist[nozero_index]
    bin_edges = bin_edges[nozero_index]
    # 计算累计频率
    cdf = np.cumsum(hist) / (msavi.size - zero_count)

    # 计算对应断点
    diff = abs(cdf - ratio * 1.0 / 100)
    laiarr_min = bin_edges[np.where(diff == diff.min())]
    diff = abs(cdf - (1 - ratio * 1.0 / 100))
    laiarr_max = bin_edges[np.where(diff == diff.min())]

    return laiarr_min, laiarr_max,msavi.mean(),msavi.std()

def findsnow(shp,tif):
    """
    """
    # 矢量栅格化参考tif
    #data = gdal.Open(tif,gdal.GA_ReadOnly)
    data = gdal.Open(tif)
    x_res = data.RasterXSize
    y_res = data.RasterYSize

    # 根据类别进行矢量转栅格处理 snow = 255; road=100; GH=25
    vector = ogr.Open(shp)
    layer = vector.GetLayer()

    # 栅格矢量化放到内存中
    targetDataSet = gdal.GetDriverByName("MEM").Create("", x_res, y_res, 1, gdal.GDT_Byte)
    targetDataSet.SetGeoTransform(data.GetGeoTransform())
    targetDataSet.SetProjection(data.GetProjection())
    band = targetDataSet.GetRasterBand(1)
    NoData_value = -999
    band.SetNoDataValue(NoData_value)
    band.FlushCache()
    gdal.RasterizeLayer(targetDataSet, [1], layer, options=["ATTRIBUTE=ND"])
    mask = band.ReadAsArray()

    # band = None
    # targetDataSet = None

    # indexsnow = np.where(mask==255)
    # indexroad = np.where(mask==100)
    # indexother = np.where(mask==50)
    # indexdapen = np.where(mask==25)

    # mask 值对应的是矢量字段中类别填充的DN值
    indexsnow = np.where(mask == 255)
    indexdapen = np.where(mask == 0)
    # indexroad = np.where(mask == 100)
    # indexother = np.where(mask == 50)

    # 构建雪体指数
    im_g =data.GetRasterBand(3).ReadAsArray()
    im_swir1 = data.GetRasterBand(11).ReadAsArray()
    nodnindex = np.where(im_g==0)
    snowvi = (im_g - im_swir1) / ((im_g + im_swir1) + 0.000001)

    # 雪、路、其他
    snow = snowvi[indexsnow]
    dapeng = snowvi[indexdapen]
    # road = snowvi[indexroad]
    # other = snowvi[indexother]

    # 设置置信度来计算数组的min\max  confidence interval
    snow_nor = (snow * 1000).astype(np.int32)
    snow_min, snowarr_max,snowarr_mean,snowarr_std = calc_confidence_interval(snow_nor, 1)
    snowmin = snow_min[0]/1000
    snowrmax = snowarr_max[0]/1000
    snowarrmean = snowarr_mean/1000
    snowarrstd  = snowarr_std/1000

    # 设置置信度来计算数组的min\max  confidence interval
    dapeng_nor = (dapeng * 1000).astype(np.int32)
    dapeng_min, dapeng_max,dapeng_mean,dapeng_std = calc_confidence_interval(dapeng_nor, 1)
    dapengmin = dapeng_min[0]/1000
    dapengmax = dapeng_max[0]/1000
    dapengmean = dapeng_mean/1000
    dapengstd = dapeng_std/1000

    # 设置置信度来计算数组的min\max  confidence interval
    # other_nor = (other * 1000).astype(np.int32)
    # other_min, other_max,other_mean,other_std = calc_confidence_interval(other_nor, 1)
    # othermin = other_min[0]/1000
    # ohtermax = other_max[0]/1000
    # othermean = other_mean/1000
    # otherstd = other_std/1000

    # 第一次阈值计算
    PRE_THDN = ((snowarrmean-0.5*snowarrstd) + (dapengmean+0.5*dapengstd))*0.5

    print("snowmin:",snowmin)
    print("PRE_THDN:",PRE_THDN)

    # 影像成像日期
    imagetimemouse = os.path.split(tif)[1].split("-")[1]
    imagetimeday = os.path.split(tif)[1].split("-")[2].split("_")[0]

    # 如果是3月份的影像，选择较大的阈值
    if imagetimemouse == "3" or imagetimemouse == "03":
        PRE_THDN = max([PRE_THDN,snowmin])
        print("3月份选择较大的阈值，new_PRE_THDN",PRE_THDN)

    # 2月中旬之前的影像，选择较小阈值
    elif imagetimemouse == "02" or imagetimemouse == "02":
        if int(imagetimeday) <= 15 :
            PRE_THDN = min([PRE_THDN, snowmin])
            print("2月中旬之前的影像，选择较小阈值new_PRE_THDN：",PRE_THDN)

    # 其他时间按正常阈值
    else:
        PRE_THDN = PRE_THDN

    snowindex = np.where(snowvi >= PRE_THDN)
    #nodataindex = np.where(snowvi == nodata)
    # nosnowindex = np.where(im_data<0.793827)
    nosnowindex = np.where((snowvi < PRE_THDN))

    snowvi[snowindex] = 200
    snowvi[nosnowindex] = 100
    snowvi[nodnindex] = 0

    tiff_driver = gdal.GetDriverByName('GTiff')
    inpath = os.path.split(tif)[0]
    imname = os.path.split(tif)[1].replace(".tif", "Thresholdpre.tif")
    out_tif = os.path.join(inpath, imname)
    out_ds = tiff_driver.Create(out_tif, x_res, y_res, 1, gdal.GDT_Float32)
    out_ds.SetProjection(data.GetProjection())
    out_ds.SetGeoTransform(data.GetGeoTransform())
    out_band = out_ds.GetRasterBand(1)
    out_band.SetNoDataValue(0)
    out_band.WriteArray(snowvi)
    out_band.FlushCache()
    return out_tif

def main(imagedir,shpdir,fidname):

    code = 200
    nodata = 0
    # 创建临时目录
    temfold = os.path.join(os.path.split(shpdir)[0], "tem")
    if os.path.isdir(temfold) == False:
        os.makedirs(temfold)
    gdal.SetConfigOption('SHAPE_ENCODING', 'GBK')
    print("start to read raster")
    dataset = gdal.Open(imagedir)
    im_project = dataset.GetProjection()
    im_tranform = dataset.GetGeoTransform()
    xsize = dataset.RasterXSize
    ysize = dataset.RasterYSize
    imdata = dataset.ReadAsArray(0, 0, xsize, ysize)
    #将元数据组合
    meta1 = [xsize, ysize, im_project, im_tranform, imdata]
    datasource = ogr.Open(shpdir,1)
    layer = datasource.GetLayer(0)
    feature_number = layer.GetFeatureCount()
    progress(0.25)
    star_end=[0,feature_number]
    out = imdata
    filldata(nodata,out, code, layer, meta1,temfold, star_end,fidname)
    progress(0.85)

if __name__ == '__main__':

    start = time.time()
    print(start)
    tif = r"E:\三江\浓江前进859\2022-02-18_4071_S2DL.tif"
    shp =r"E:\三江\浓江前进859\浓江前进859\sampleNJ3.shp"
    out_tif = findsnow(shp, tif) # 输出的是二值化后的影像，便于判断清雪区域是否正确
    fidname = "0307QXb"
    main(out_tif,shp,fidname)
    end = time.time()
    print("耗时：",end -start)
