#############################################################################
##
#W  mapper.gd                  GAP library                   Alexander Hulpke
##
##
#Y  Copyright (C) 2016 The GAP Group
##
##  This file declares the operations for mappers, objects that can be used
##  as underlying machinery for homomorphisms but do not aim to be algebraic
##  objects themselves.
##

# the domain of a mapper might be only defined implicitly in the implementation.
DeclareCategory("IsMapper",IsObject);

# a mapper that uses group properties
DeclareCategory("IsGroupMapper",IsMapper);

# all mappers are in the same family (as we do not declare domain/range)
BindGlobal("MapperFamily",NewFamily( "MapperFamily", IsMapper ));

#############################################################################
##
#O  ImagesRepresentative(<mapr>,<elm>) . one image of elm under a mapper
##
##  <#GAPDoc Label="ImagesRepresentative">
##  <ManSection>
##  <Oper Name="ImagesRepresentative" Arg='mapr,elm'/>
##
##  <Description>
##  returns an image when applying a mapper to an element of its (implicitly
##  defined) domain. If <A>elm</A> is not it this domain the behavior is
##  undefined.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareOperation( "ImagesRepresentative", [ IsMapper, IsObject ] );


#############################################################################
##
#A  InverseMapper( <Mapr> )
##
##  <#GAPDoc Label="InverseMapper">
##  <ManSection>
##  <Attr Name="InverseMapper" Arg='Mapr'/>
##
##  <Description>
##  The inverse of a mapper <A>M</A>is a mapper <A>I</A> that satisfies 
##  that <A>M*I*M</A> gives the same image as <A>M</A>. This implies that
##  inverse mappers need to be applicable at least to the range of original
##  mapper, and that inverse mappers may store each other mutually.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareAttribute( "InverseMapper", IsMapper );

#############################################################################
##
#A  FamilySourceMapper( <Mapr> )
##
##  <#GAPDoc Label="FamilySourceMapper">
##  <ManSection>
##  <Attr Name="FamilySourceMapper" Arg='Mapr'/>
##
##  <Description>
##  This attribute, if set, holds the family of objects to which the mapper
##  is applicable. It can be used as first step in the membership test.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareAttribute( "FamilySource", IsMapper );

#############################################################################
##
#P  CanTestApplicabilityMapper( <Mapr> )
##
##  <#GAPDoc Label="CanTestApplicabilityMapper">
##  <ManSection>
##  <Prop Name="CanTestApplicabilityMapper" Arg='Mapr'/>
##
##  <Description>
##  A mapper for which this property has been set can be applied to
##  arbitrary elements of its <K>FamilySourceMapper</K> and will either return
##  a valid image, or <K>fail</K>, if the element is not in the implicitly
##  defined domain.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareProperty( "CanTestApplicabilityMapper", IsMapper );

#############################################################################
##
#A  MappingGeneratorsImages(<mapr>)
#O  MappingGeneratorsImages(<mapr>,<dom>)
##
##  <#GAPDoc Label="MappingGeneratorsImages">
##  <ManSection>
##  <Attr Name="MappingGeneratorsImages" Arg='mapr,dom'/>
##  <Oper Name="MappingGeneratorsImages" Arg='mapr,dom'/>
#
##  <Description>
##  For a mapper and a domain, this operation returns a list of generators
##  for the domain, and a corresponding list of images under the mapper.
##  In the attribute version, generators and images are stored to define the
##  mapper.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareOperation( "MappingGeneratorsImages", [IsMapper,IsObject] );
DeclareAttribute( "MappingGeneratorsImages", IsMapper );

#############################################################################
##
#A  MultiplicativeKernelMapper(<mapr>,<dom>)
##
##  <#GAPDoc Label="MultiplicativeKernelMapper">
##  <ManSection>
##  <Oper Name="MultiplicativeKernelMapper" Arg='mapr,dom'/>
#
##  <Description>
##  For a mapper preserving multiplication and a domain with an
##  associative multiplication, this operation returns the kernel of the
##  homomorphism defined by the mapper on this domain.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareOperation( "MultiplicativeKernelMapper", [IsMapper,IsObject] );

#############################################################################
##
#A  IsSingleValuedMapper(<mapr>,<dom>)
##
##  <#GAPDoc Label="IsSingleValuedMapper">
##  <ManSection>
##  <Oper Name="IsSingleValuedMapper" Arg='mapr,dom'/>
#
##  <Description>
##  Tests whether the mapper, restricted to the given domain, is single
##  valued.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareOperation( "IsSingleValuedMapper", [IsMapper,IsObject] );

# --- create mappers

#############################################################################
##
#O  GroupMapperByImages(<G>,<gens>,<imgs>)  . . . create mapper
##
##  <#GAPDoc Label="GroupMapperByImages">
##  <ManSection>
##  <Oper Name="GroupMapperByImages" Arg='G, gens, imgs'/>
##
##  <Description>
##  Given a group <A>G</A> and a set of generators <A>gens</A> for this
##  group, the operation Constructs a mapper that maps <A>gens</A> to
##  <A>imgs</A>.
##  If the span of <A>gens</A> differs from <A>G</A>, the result is
##  undefined.
##  Group and generators are given to allow for different choice of
##  generators than stored, but use of group information.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareOperation( "GroupMapperByImages",
    [ IsGroup, IsList, IsList ] );

#############################################################################
##
#O  GroupMapperByAction(<dom>,<act>)  . . . create mapper
##
##  <#GAPDoc Label="GroupMapperByAction">
##  <ManSection>
##  <Oper Name="GroupMapperByAction" Arg='dom, act'/>
##
##  <Description>
##  Creates a mapper that will map elements to the permutation
##  of <A>dom</A> (interpreted as a sorted list) obtained by acting with the
##  function <A>act</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareOperation( "GroupMapperByAction",
    [  IsCollection, IsFunction ] );

#############################################################################
##
#O  GroupMapperByFunction(<fct>)  . . . create mapper
##
##  <#GAPDoc Label="GroupMapperByFunction">
##  <ManSection>
##  <Oper Name="GroupMapperByFunction" Arg='fct'/>
##
##  <Description>
##  Creates a mapper that will map elements by applying the function
##  <A>fct</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareOperation( "GroupMapperByFunction", [ IsFunction ] );


# --- mapper composition

#############################################################################
##
#O  CompositionMapper2(<map2>,<map1>)  . . . composition of mappers
##
##  <#GAPDoc Label="CompositionMapper2">
##  <ManSection>
##  <Oper Name="CompositionMapper2" Arg='map2, map1'/>
##
##  <Description>
##  <Ref Func="CompositionMapper2"/> returns a mapper that applies as a
##  sequence of first applying <A>map1</A> and then applying <A>map2</A>.
##  (Note the reverse ordering of arguments in the composition via
##  the multiplication <Ref Func="\*"/>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareOperation( "CompositionMapper2",
    [ IsMapper, IsMapper ] );


#############################################################################
##
#R  IsCompositionMapperRep( <map> )
##
##  <#GAPDoc Label="IsCompositionMapperRep">
##  <ManSection>
##  <Filt Name="IsCompositionMapperRep" Arg='map' Type='Representation'/>
##
##  <Description>
##  Mappers in this representation are stored as composition of two
##  mappers.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##
DeclareRepresentation( "IsCompositionMappingRep",
    IsMapper and IsAttributeStoringRep, [ "map1", "map2" ] );

# --- Representations for particular mappers

DeclareRepresentation( "IsMapperRep", IsMapper and IsAttributeStoringRep,[]);

DeclareRepresentation( "IsPermGroupMapper", IsMapperRep and IsGroupMapper,[]);
DeclareRepresentation( "IsToPermGroupMapper", IsMapperRep and IsGroupMapper,[]);

