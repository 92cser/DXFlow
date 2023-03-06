# 关于DXFlow

各位老师、同学，大家好！我是王东旭，2022年毕业于中国海洋大学，现就职于宁波大学海洋工程研究院。DXFlow是由我开发的基于OpenFOAM的数值波浪水槽程序，与waves2foam和OLAFlow类似，但也有一些自己的风格，如改进的源造波方法、海绵层消波方法、重叠网格和isoAdvector等，大多都有论文与之对应（关于我发表的论文，可在[这里](http://ioe.nbu.edu.cn/info/1040/1124.htm)查看）。经过深思熟虑，我决定先将beta版开源，供大家交流学习。

当前版本的DXFlow需基于OFv2006编译，没有做多版本的适应性调整。同时，我认为程序的功能还不够丰富，故暂时没有提供用户手册（没什么内容，大家用用就都会了）。

我计划于每6-8个月更新一次，release新的版本。目前，正在开发的方向有：

1. Ghost Fluid method (基本完成，需雕琢细节)
2. New mesh intepolation scheme of overset grid（完成，通过计算评估结果）
3. velocity input method for wave / current generaiton （不是新算法，有时间想起来了就搞一下）

欢迎各位老师、同学试用DXFlow，并基于其开展工作！如果能与我开展合作就更好了：）

最后，科研的关键在于交流，故步自封不能长久，闭门造车难以发展，希望DXFlow能对国内开源事业的发展起到推动作用。

请注意：DXFlow只供学习交流，同时由于本人水平有限，程序一定也许多不足。因此，

## 请勿将DXFlow用于商业用途！请勿将DXFlow用于商业用途！请勿将DXFlow用于商业用途！


（先在国内开源吧，国外友人想用就用，下面是个极简版的introduction）
# Description
DXFlow is an open-sourced CFD toolbox based on OpenFOAM 2006. It is like [waves2Foam](https://openfoamwiki.net/index.php/Contrib/waves2Foam) and [OLAFlow](https://openfoamwiki.net/index.php/Contrib/olaFlow) ([IHFoam](https://openfoamwiki.net/index.php/Contrib/IHFOAM)) but with the following interest features:

* (Line-shaped) mass source wavemaker with velocity adjustment inside source region;

* First- and second-order piston type wavemaker (with active absorption);

* Combination of overset grid method and isoAdvector.

# DXFlow download and compilation

* download and compile OpenFOAM 2006 at www.openfoam.com

* ```$ git clone https://github.com/92cser/DXFlow.git```

* ```$ ./Allwmake```

# Note
DXFlow is totally free and only for study use. See the GNU General Public License for more details.
