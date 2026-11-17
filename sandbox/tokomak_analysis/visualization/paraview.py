"""Run with pvpython: import the multiblock, colour and optionally cut away."""
import argparse


def main():
    parser=argparse.ArgumentParser();parser.add_argument('vtm');parser.add_argument('--cutaway',action='store_true');parser.add_argument('--screenshot',default='tokamak.png')
    parser.add_argument('--color',choices=['material','component'],default='material')
    args=parser.parse_args()
    from paraview.simple import XMLMultiBlockDataReader,Show,ColorBy,GetActiveViewOrCreate,Clip,ResetCamera,SaveScreenshot
    reactor=XMLMultiBlockDataReader(FileName=[args.vtm])
    view=GetActiveViewOrCreate('RenderView');display=Show(reactor,view);ColorBy(display,('CELLS',args.color+'_id'))
    if args.color=='material':
        import json
        from pathlib import Path
        from paraview.simple import GetColorTransferFunction
        colors=json.loads((Path(args.vtm).parent/'material_colors.json').read_text())
        lut=GetColorTransferFunction('material_id')
        lut.RGBPoints=[v for item in colors for v in [item['material_id'],*[x/255 for x in item['rgb']]]]
    if args.cutaway:
        display.Visibility=0;clip=Clip(Input=reactor);clip.ClipType.Normal=[0,1,0];Show(clip,view)
    ResetCamera(view);SaveScreenshot(args.screenshot,view,ImageResolution=[1600,1200])


if __name__=='__main__':main()
