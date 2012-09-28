<? include_once("BasicTypes.inc.php");

class FrogFEM extends SOFAObject{
public function FrogFEM(){
		$ab=new AABB(new Coord(-10,-3,-7),new Coord(8,2.5,7));
		$this->setAABB($ab);
}
public function printObj(){
?>
	<Node name="FrogFEM" processor="0">
		<EulerImplicitSolver rayleighStiffness="0.1" rayleighMass="0" />
    <ParallelCGLinearSolver iterations="25" threshold="0.000000000001" tolerance="0.000001" />
		<SparseGrid n="6 6 6" filename="mesh/frog_body.obj"/>
		<MechanicalObject dx="<?=$this->dx?>" dy="<?=$this->dy?>" dz="<?=$this->dz?>" rx="<?=$this->rx?>" rz="<?=$this->rz?>" ry="<?=$this->ry?>" />
		<UniformMass totalmass="5" />
		<HexahedronFEMForceField name="FEM" youngModulus="5000" poissonRatio="0.3"  method="polar"/>
		<Node name="Visu1">
			<OglModel name="VisualBody" fileMesh="mesh/frog_body.obj" normals="0" dx="<?=$this->dx?>" dy="<?=$this->dy?>" dz="<?=$this->dz?>" rx="<?=$this->rx?>" rz="<?=$this->rz?>" ry="<?=$this->ry?>"  color="0.17 0.70 0.05" />
			<BarycentricMapping object1="../.." object2="VisualBody" />
		</Node>
		<Node name="Visu2">
			<OglModel name="VisualEyes" fileMesh="mesh/frog_eyes.obj" normals="0" dx="<?=$this->dx?>" dy="<?=$this->dy?>" dz="<?=$this->dz?>" rx="<?=$this->rx?>" rz="<?=$this->rz?>" ry="<?=$this->ry?>" color="0.04 0.19 0.52" />
			<BarycentricMapping object1="../.." object2="VisualEyes" />
		</Node>
		<Node name="Visu3">
			<OglModel name="VisualEyebrows" fileMesh="mesh/frog_eyebrows.obj" normals="0"  dx="<?=$this->dx?>" dy="<?=$this->dy?>" dz="<?=$this->dz?>" rx="<?=$this->rx?>" rz="<?=$this->rz?>" ry="<?=$this->ry?>"  color="0.44 0.43 0.00" />
			<BarycentricMapping object1="../.." object2="VisualEyebrows" />
		</Node>
		<Node name="Visu4">
			<OglModel name="VisualLips" fileMesh="mesh/frog_lips.obj" normals="0"  dx="<?=$this->dx?>" dy="<?=$this->dy?>" dz="<?=$this->dz?>" rx="<?=$this->rx?>" rz="<?=$this->rz?>" ry="<?=$this->ry?>"  color="0.47 0.25 0.03" />
			<BarycentricMapping object1="../.." object2="VisualLips" />
		</Node>
		<Node name="CollisionNode">
	  	<MeshLoader name="collisionLoader"  filename="mesh/frog-push25.obj"/>
			<Mesh />
			<MechanicalObject  dx="<?=$this->dx?>" dy="<?=$this->dy?>" dz="<?=$this->dz?>" rx="<?=$this->rx?>" rz="<?=$this->rz?>" ry="<?=$this->ry?>" />
			<Triangle />
			<BarycentricMapping />
		</Node>
	</Node>
<?}
};
?>
