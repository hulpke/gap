#############################################################################
##
#W  mapper.gd                  GAP library                   Alexander Hulpke
#W                                                 Ákos Seress, Heiko Theißen
##
##
#Y  Copyright (C) 2016 The GAP Group
##
##  This file contains the implementation for mappers for permutation groups.
##

InstallMethod(GroupMapperByImages,"perm group",IsFamFamX,
  [IsPermGroup,IsList,IsList],0,
function(G,gens,imgs)
local filter,obj;
  if Length(imgs)=0 then
    TryNextMethod(); # if empty images need to give one
  fi;
  Assert(1,Length(gens)=Length(imgs));
  Assert(2,ForAll(gens,x->x in G));
  Assert(2,G=Group(gens));
  filter:=IsPermGroupMapper;
  if ForAll(imgs,IsPerm) then
    filter:=filter and IsToPermGroupMapper;
  fi;
  obj:=Objectify(NewType(MapperFamily,filter),rec(group:=G,oneimg:=One(imgs[1])));
  SetMappingGeneratorsImages(obj,[gens,imgs]);
  return obj;
end);

InstallMethod(ViewObj,"mapper with gensimgs",true,
  [IsMapper and HasMappingGeneratorsImages],0,
function(m)
  m:=MappingGeneratorsImages(m);
  ViewObj(m[1]);
  PrintObj("=>");
  ViewObj(m[2]);
end);


#############################################################################
##
#M  StabChainMutable( <hom> ) . . . . . . . . . . . . . . for perm mappers
##
InstallOtherMethod( StabChainMutable, "perm mapper",  true,
        [ IsPermGroupMapper ], 0,
    function( map )
    local   S,
            rnd,        # list of random elements of '<hom>.source'
            rne,        # list of the images of the elements in <rnd>
            rni,        # index of the next random element to consider
            elm,        # one element in '<hom>.source'
            img,        # its image
            size,       # size of the stabilizer chain constructed so far
            stb,        # stabilizer in '<hom>.source'
            bpt,        # base point
            two,        # power of two
            trivgens,   # trivial generators and their images, must be
            trivimgs,   #   entered into every level of the chain
            mapi,
            i, T,  # loop variables
            orb,
            orbf,       # indicates with which generator the image was obtained
            dict,
            cnt,
            short,
            FillTransversalShort,
            BuildOrb,
            AddToStbO,
            maxstor,
            gsize,
            l;  # position

    # Add to short word orbit fct.
    AddToStbO:=function(o,dict,e,w)
    local i;
      #Print("add length ",Length(UnderlyingElement(w)),"\n");
      i:=LookupDictionary(dict,e);
      if i<>fail then
        if Length(o[i][2])>Length(w) then
          o[i]:=Immutable([e,w]);
          return 0;
        fi;
        return 1;
      else
        Add(o,Immutable([e,w]));
        AddDictionary(dict,e,Length(o));
        return 0;
      fi;

#      if l<>Fail then
#      for i in [1..Length(o)] do
#       if o[i][1]=e then
#         if Length(o[i][2])>Length(w) then
#           o[i]:=Immutable([e,w]);
#         fi;
#         return;
#       fi;
#      od;
#      Add(o,Immutable([e,w]));
    end;

    # build short words by an orbit algorithm on genimg
    BuildOrb:=function(genimg)
    local a,orb,dict,orbf,T,elm,img,i,n;
      if Length(genimg[1])>0 then
	a:=genimg[1][1];
      else
	a:=One(map!.group);
      fi;
      dict:=NewDictionary(a,false);
      a:=One(map!.group);
      AddDictionary(dict,a);
      orb:=[Immutable([a,map!.oneimg])];
      orbf:=[0];
      i:=1;
      n:=Length(genimg[1]);
      while Length(orb)<maxstor and i<=Length(orb) do
        for T in [1..n] do
          if orbf[i]<>-T then
            elm:=orb[i][1]*genimg[1][T];
            if not KnowsDictionary(dict,elm) then
              # new rep found
              img:=orb[i][2]*genimg[2][T];
              AddDictionary(dict,elm);
              Add(orb,Immutable([elm,img]));
              Add(orbf,T);
            fi;
          fi;
          if orbf[i]<>T then
            elm:=orb[i][1]/genimg[1][T];
            if not KnowsDictionary(dict,elm) then
              # new rep found
              img:=orb[i][2]/genimg[2][T];
              AddDictionary(dict,elm);
              Add(orb,Immutable([elm,img]));
              Add(orbf,-T);
            fi;
          fi;
        od;
        i:=i+1;
      od;
      return orb;
    end;

    mapi:=MappingGeneratorsImages(map);

    # do products build up? (Must we prefer short words?)
    short:=(IsAssocWordCollection(mapi[2]) 
         or IsElementOfFpGroupCollection(mapi[2]))
            and ValueOption("noshort")<>true;

    if short then
      # compute how many perms we permit to store?
      maxstor:=LargestMovedPoint(mapi[1])+1;
      if maxstor>65535 then
        maxstor:=maxstor*2; # perms need twice as much memory
      fi;
      maxstor:=Int(40*1024^2/maxstor); # allocate at most 40MB to the perms
      # but don't be crazy 
      maxstor:=Minimum(maxstor,
                 Size(map!.group)/10,
                 500*LogInt(Size(map!.group),2),
                 25000); 

      # fill transversal with elements that are short words
      # This is similar to Minkwitz' approach and produces much shorter
      # words when decoding.
      FillTransversalShort:=function(stb,size)
      local l,i,bpt,m,elm,wrd,z,j,dict,fc,mfc;
        mfc:=Minimum(maxstor*10,gsize/size);
        bpt:=stb.orbit[1];
        stb.norbit:=ShallowCopy(stb.orbit);
        # fill transversal with short words
        for l in stb.orb do
          i:=bpt/l[1];
          if not i in stb.norbit then
            Add(stb.norbit,i);
            stb.transversal[i]:=l[1];
            stb.transimages[i]:=l[2];
          fi;
          i:=bpt^l[1];
          if not i in stb.norbit then
            Add(stb.norbit,i);
            stb.transversal[i]:=Inverse(l[1]);
            stb.transimages[i]:=Inverse(l[2]);
          fi;
        od;
        stb.stabilizer.orb:=Filtered(stb.orb,i->bpt^i[1]=bpt);
        dict:=NewDictionary(stb.stabilizer.orb[1][1],true);
        for l in [1..Length(stb.stabilizer.orb)] do
          AddDictionary(dict,stb.stabilizer.orb[l][1],l);
        od;
        l:=1;
        fc:=1;
        maxstor:=Minimum(maxstor,QuoInt(5*gsize,size));
        if maxstor<1000 then 
          maxstor:=Maximum(maxstor,Minimum(QuoInt(gsize,size),1000));
        fi;
        #Print(maxstor," ",gsize/size,"<\n");
        while Length(stb.stabilizer.orb)*5<maxstor and l<=Length(stb.orb)
          and fc<mfc do
          # add schreier gens
          elm:=stb.orb[l][1];
          wrd:=stb.orb[l][2];
          for z in [1,2] do
            if z=2 then
              elm:=elm^-1;
              wrd:=wrd^-1;
            fi;
            i:=bpt^elm;
            for j in stb.orb do
              if bpt^j[1]=i then
                fc:=fc+AddToStbO(stb.stabilizer.orb,dict,elm/j[1],wrd/j[2]);
              elif i^j[1]=bpt then
                fc:=fc+AddToStbO(stb.stabilizer.orb,dict,elm*j[1],wrd*j[2]);
              fi;
            od;
          od;
          l:=l+1;
        od;

        Unbind(stb.orb);
        Unbind(stb.norbit);
        stb:=stb.stabilizer;
        #Print("|o|=",Length(stb.orb),"\n");
        # is there too little left? If yes, extend!
        if Length(stb.orb)*20<maxstor then
          stb.orb:=BuildOrb([List(stb.orb,i->i[1]),
                             List(stb.orb,i->i[2])]);
        fi;
#Print(bpt,":",Length(stb.orb),"\n");
      end;
    else
      FillTransversalShort:=Ignore;
    fi;

    # initialize the random generators
    two := 16;
    rnd := ShallowCopy( mapi[1] );
    for i  in [Length(rnd)..two]  do
        Add( rnd, One( map!.group ) );
    od;
    rne := ShallowCopy( mapi[2] );
    for i  in [Length(rne)..two]  do
        Add( rne, map!.oneimg );
    od;
    rni := 1;

    S := EmptyStabChain( [  ], One( map!.group ),
                         [  ], map!.oneimg );
    if short then
      S.orb:=BuildOrb(mapi);
    fi;

    # initialize the top level
    bpt:=fail;
    if short then
      bpt:=DoShortwordBasepoint(S.orb);
    fi;
    if bpt=fail then;
      bpt := SmallestMovedPoint( mapi[1] );
      if bpt = infinity  then
          bpt := 1;
      fi;
    fi;
    InsertTrivialStabilizer( S, bpt );
    # the short words usable on this level
    gsize:=Size(map!.group);
    FillTransversalShort(S,1);
    
    # Extend  orbit and transversal. Store  images of the  identity for other
    # levels.
    AddGeneratorsGenimagesExtendSchreierTree( S, mapi[1], mapi[2] );
    trivgens := [  ];  trivimgs := [  ];
    for i  in [ 1 .. Length( mapi[1] ) ]  do
        if mapi[1][ i ] = One( map!.group )  then
            Add( trivgens, mapi[1][ i ] );
            Add( trivimgs, mapi[2][ i ] );
        fi;
    od;

    # get the size of the stabilizer chain
    size := Length( S.orbit );

    # create new elements until we have reached the size
    while size <> gsize  do

        # try random elements
        elm := rnd[rni];
        img := rne[rni];
        i := Random( [ 1 .. Length( mapi[1] ) ] );
        rnd[rni] := rnd[rni] * mapi[1][i];
        rne[rni] := rne[rni] * mapi[2][i];
        rni := rni mod two + 1;

        # divide the element through the stabilizer chain
        stb := S;
        bpt := BasePoint( stb );
        while     bpt <> false
              and elm <> stb.identity
              and Length( stb.genlabels ) <> 0  do
            i := bpt ^ elm;
            if IsBound( stb.translabels[ i ] )  then
                while i <> bpt  do
                    img := img * stb.transimages[ i ];
                    elm := elm * stb.transversal[ i ];
                    i := bpt ^ elm;
                od;
                stb := stb.stabilizer;
                bpt := BasePoint( stb );
            else
                bpt := false;
            fi;
        od;

        # if the element was not in the stabilizer chain
        if elm <> stb.identity  then

          # if this stabilizer is trivial add an new level
          if not IsBound( stb.stabilizer )  then
            l:=fail;
            if short and IsBound(stb.orb) then
              l:=DoShortwordBasepoint(stb.orb);
            fi;
            if l=fail then
              l:=SmallestMovedPoint(elm);
            fi;
            InsertTrivialStabilizer( stb, l );
            AddGeneratorsGenimagesExtendSchreierTree( stb,
                    trivgens, trivimgs );
            # the short words usable on this level
            FillTransversalShort(stb,size);
          fi;

#         if short then
#           l:=LookupDictionary(dict,elm);
#           if l<>fail then
#             img:=l;
#           fi;
#         fi;

          # extend the Schreier trees above level `stb'
          T := S;
          repeat
            T := T.stabilizer;
            size := size / Length( T.orbit );
            AddGeneratorsGenimagesExtendSchreierTree( T, [elm], [img] );
            size := size * Length( T.orbit );
          until T.orbit[ 1 ] = stb.orbit[ 1 ];

        fi;

    od;
    
    return S;
end );

InstallOtherMethod( StabChainMutable, "perm to perm mapper",  true,
        [ IsPermGroupMapper and IsToPermGroupMapper], 0,
    function( map )
    local   options,    # options record for stabilizer construction
            n,  
            k,
            i,
            a,b,
            longgens,
            longgroup,
            conperm,
            conperminv,
            mapi,
            op;
    
    mapi:=MappingGeneratorsImages(map);
    n := LargestMovedPoint( mapi[1] );
    k := LargestMovedPoint( mapi[2] );
    
    # collect info for options
    options := rec();
    
    # random or deterministic
    if   IsBound( StabChainOptions( Parent( Source( hom ) ) ).random )  then
        options.randomSource :=
          StabChainOptions( Parent( map!.group ) ).random;
    elif IsBound( StabChainOptions( map!.group ).random )  then
        options.randomSource := StabChainOptions( map!.group ).random;
    #elif IsBound( StabChainOptions( PreImagesRange( hom ) ).random )  then
    #    options.randomSource := StabChainOptions( PreImagesRange( hom ) ).random;
    else
        options.randomSource := DefaultStabChainOptions.random;
    fi;

    if   IsBound(map!.range) and IsBound( StabChainOptions( Parent( map!.range ) ).random )  then
        options.randomRange :=
          StabChainOptions( Parent( map!.range ) ).random;
    elif IsBound(map!.range) and IsBound( StabChainOptions( map!.range ).random )  then
        options.randomRange := StabChainOptions( map!.range ).random;
    else
      options.randomRange := DefaultStabChainOptions.random;
    fi;
    options.random := Minimum(options.randomSource,options.randomRange);
    Error("TODO");

    # if IsMapping, try to extract info from source
    if Tester( IsMapping )( hom )  and  IsMapping( hom )  then
        if   HasSize( Source( hom ) )  then
            options.size := Size( Source( hom ) );
        elif HasSize( PreImagesRange( hom ) )  then
            options.size := Size( PreImagesRange( hom ) );
        fi;
        if not IsBound( options.size )
           and HasSize( Parent( Source( hom ) ) )  then
            options.limit := Size( Parent( Source( hom ) ) );
        fi;
        if   IsBound( StabChainOptions( Source( hom ) ).knownBase )  then
            options.knownBase := StabChainOptions( Source( hom ) ).knownBase;
        elif IsBound( StabChainOptions( PreImagesRange( hom ) ).knownBase )
          then
            options.knownBase := StabChainOptions( PreImagesRange( hom ) ).
                                 knownBase;
        elif HasBaseOfGroup( Source( hom ) )  then
            options.knownBase := BaseOfGroup( Source( hom ) );
        elif HasBaseOfGroup( PreImagesRange( hom ) )  then
            options.knownBase := BaseOfGroup( PreImagesRange( hom ) );
        elif IsBound( StabChainOptions( Parent( Source( hom ) ) ).knownBase )
          then
            options.knownBase :=
              StabChainOptions( Parent( Source( hom ) ) ).knownBase;
        elif HasBaseOfGroup( Parent( Source( hom ) ) )  then
            options.knownBase := BaseOfGroup( Parent( Source( hom ) ) );
        fi;

    # if not IsMapping, settle for less
    else
        if   HasSize( Source( hom ) )  then
            options.limitSource := Size( Source( hom ) );
        elif HasSize( PreImagesRange( hom ) )  then
            options.limitSource := Size( PreImagesRange( hom ) );
        elif HasSize( Parent( Source( hom ) ) )  then
            options.limitSource := Size( Parent( Source( hom ) ) );
        fi;
        if   IsBound( StabChainOptions( Source( hom ) ).knownBase )  then
            options.knownBaseSource :=
              StabChainOptions( Source( hom ) ).knownBase;
        elif IsBound( StabChainOptions( PreImagesRange( hom ) ).knownBase )
          then
            options.knownBaseSource :=
              StabChainOptions( PreImagesRange( hom ) ).knownBase;
        elif IsBound( StabChainOptions( Parent( Source( hom ) ) ).knownBase )
          then
            options.knownBaseSource :=
                StabChainOptions( Parent( Source( hom ) ) ).knownBase;
        fi;

        # if we have info about source, try to collect info about range
        if IsBound( options.limitSource ) then 
            if   HasSize( Range( hom ) )  then
                options.limitRange := Size( Range( hom ) );
            elif HasImagesSource(hom) and HasSize( ImagesSource( hom ) )  then
                options.limitRange := Size( ImagesSource( hom ) );
            elif HasSize( Parent( Range( hom ) ) )  then
                options.limitRange := Size( Parent( Range( hom ) ) );
            fi;
            if IsBound( options.limitRange ) then 
                options.limit := options.limitSource * options.limitRange;
            fi;
        fi;
        if IsBound( options.knownBaseRange ) then 
            if   IsBound( StabChainOptions( Range( hom ) ).knownBase )  then
                options.knownBaseRange :=
                  StabChainOptions( Range( hom ) ).knownBase;
            elif IsBound( StabChainOptions( PreImagesRange( hom ) ).
                    knownBase )  then
                options.knownBaseRange :=
                  StabChainOptions( PreImagesRange( hom ) ).knownBase;
            elif IsBound( StabChainOptions( Parent( Range( hom ) ) )
                    .knownBase )
              then
                options.knownBaseRange :=
                    StabChainOptions( Parent( Range( hom ) ) ).knownBase;
            fi;
            if IsBound( options.knownBaseRange ) then 
                options.knownBase := Union( options.knownBaseSource,
                                            options.knownBaseRange + n );
            fi;
        fi;

    fi; # if IsMapping

    options.base := [1..n];

    # create concatenation of perms in hom.generators, hom.genimages
    longgens := [];
    conperm := MappingPermListList([1..k],[n+1..n+k]);
    conperminv := conperm^(-1);
    for i in [1..Length(mapi[1])] do
      # this is necessary to remove spurious points if the permutations are
      # not internal
      a:=mapi[1][i];
      b:=mapi[2][i];
      if not IsInternalRep(a) then
        a:=RestrictedPermNC(a,[1..n]);
      fi;
      if not IsInternalRep(b) then
        b:=RestrictedPermNC(b,[1..k]);
      fi;
      longgens[i] := a * (b ^ conperm); 
    od;
    longgroup :=  GroupByGenerators( longgens, One( map!.group ) );
    Error("TODO");
    for op  in [ PreImagesRange, ImagesSource ]  do
        if      Tester(op)(hom) and HasIsSolvableGroup( op( hom ) )
           and not IsSolvableGroup( op( hom ) )  then
            SetIsSolvableGroup( longgroup, false );
            break;
        fi;
    od;

    MakeStabChainLong( hom, StabChainOp( longgroup, options ),
           [ 1 .. n ], One( map!.group ), conperminv, hom,
           CoKernelOfMultiplicativeGeneralMapping );
    
    if  NrMovedPoints(longgroup)<=10000 and
       (not HasInverseGeneralMapping( hom )
       or not HasStabChainMutable( InverseGeneralMapping( hom ) )
       or not HasKernelOfMultiplicativeGeneralMapping( hom ) 
       )then
	Error("TODO");
        MakeStabChainLong( InverseGeneralMapping( hom ),
                StabChainOp( longgroup, [ n + 1 .. n + k ] ),
                [ n + 1 .. n + k ], conperminv, One( map!.group ), hom,
                KernelOfMultiplicativeGeneralMapping );
    fi;

    return StabChainMutable( map );
end );


#############################################################################
##
#M  ImagesRepresentative( <hom>, <elm> )  . . . . . . . . for perm group homs
##
InstallMethod( ImagesRepresentative, "perm group hom",FamSourceEqFamElm,
        [ IsPermGroupMapper,
          IsMultiplicativeElementWithInverse ],
function( hom, elm )
local   S,img,img2;
  S := StabChainMutable( hom );
  img := ImageSiftedBaseImage( S, OnTuples( BaseStabChain( S ), elm ),
		  S.idimage, OnRight );
  if IsPerm( img ) then
    # perm to perm group -- special treatment
    if IsInternalRep( img ) then
      TRIM_PERM( img, LargestMovedPoint( Range( hom ) ) );
    else
      img:=RestrictedPermNC(img,[1..LargestMovedPoint(Range(hom))]);
    fi;
  elif false and IsAssocWord(img) or IsElementOfFpGroup(img) then
    # try the inverse as well -- it might give a shorter image
    img2:= ImageSiftedBaseImage( S, List(BaseStabChain(S),i->i/elm),
		  S.idimage, OnRight );
    if Length(UnderlyingElement(img2))<Length(UnderlyingElement(img)) then
      return img2;
    fi;
  fi;
  return img^-1;
end );

